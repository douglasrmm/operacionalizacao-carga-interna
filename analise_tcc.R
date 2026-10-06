# ==============================================================================
# TCC - Sensibilidade das inferencias de carga de treinamento a operacionalizacao
# da relacao entre PSE e duracao
#
# Dados: Hu et al. (2024), PLOS ONE, DOI 10.1371/journal.pone.0288345
# Entrada: hu_rpe_analitico.csv
# Saida: pasta resultados/
# ==============================================================================

# 1. PACOTES E CONFIGURACAO ----------------------------------------------------
pacotes <- c("lme4", "lmerTest", "ordinal", "emmeans", "pbkrtest")
ausentes <- pacotes[!vapply(pacotes, requireNamespace, logical(1), quietly = TRUE)]
if (length(ausentes) > 0) stop("Instale antes de executar: ", paste(ausentes, collapse = ", "))

library(lme4)
library(lmerTest)
library(ordinal)
library(emmeans)
library(pbkrtest)
emm_options(lmer.df = "kenward-roger")

arquivo_dados <- "hu_rpe_analitico.csv"
if (!file.exists(arquivo_dados)) stop("Arquivo nao encontrado: ", arquivo_dados)
if (!dir.exists("resultados")) dir.create("resultados")

# 2. IMPORTACAO E PREPARACAO --------------------------------------------------
dados <- read.csv(arquivo_dados, stringsAsFactors = FALSE)

colunas_esperadas <- c("athlete", "position", "week", "model", "day", "rpe",
                       "duration", "srpe", "session")
if (!all(colunas_esperadas %in% names(dados))) {
  stop("O banco nao contem todas as colunas esperadas.")
}

dados$athlete  <- factor(dados$athlete)
dados$position <- factor(dados$position, levels = c(1, 2), labels = c("Forwards", "Backs"))
dados$model    <- factor(dados$model, levels = c(1, 2), labels = c("Model_I", "Model_II"))
dados$day      <- factor(dados$day, levels = c(1, 2, 3), labels = c("StD", "EnD", "SpD"))
dados$week     <- factor(dados$week)
dados$session  <- factor(dados$session)
dados$rpe_ord  <- ordered(dados$rpe)

# Auditoria minima do banco congelado
stopifnot(nrow(dados) == 780)
stopifnot(nlevels(dados$athlete) == 26)
stopifnot(nlevels(dados$week) == 10)
stopifnot(nlevels(dados$session) == 30)
stopifnot(!anyNA(dados[colunas_esperadas]))
stopifnot(all(abs(dados$srpe - dados$rpe * dados$duration) < 1e-10))
stopifnot(!any(duplicated(dados[c("athlete", "week", "day")])))
stopifnot(all(table(dados$athlete) == 30))

# Uma duracao por sessao
dur_por_sessao <- aggregate(duration ~ session, dados, function(x) length(unique(x)))
stopifnot(all(dur_por_sessao$duration == 1))

dados_sessao <- dados[!duplicated(dados$session), c("session", "week", "model", "day", "duration")]
stopifnot(nrow(dados_sessao) == 30)

# 3. REPRODUCAO DESCRITIVA DAS CELULAS ---------------------------------------
descritivos_srpe <- aggregate(
  srpe ~ position + model + day,
  dados,
  function(x) c(n = length(x), media = mean(x), dp = sd(x))
)
descritivos_srpe <- do.call(data.frame, descritivos_srpe)
names(descritivos_srpe)[4:6] <- c("n", "media", "dp")

# 4. DIAGNOSTICO CONVENCIONAL DA sRPE ----------------------------------------
modelo_srpe_lm <- lm(srpe ~ position * model * day, data = dados)
residuos_lm_pad <- rstandard(modelo_srpe_lm)
shapiro_srpe <- shapiro.test(residuos_lm_pad)

# 5. MODELO DE REFERENCIA: PSE x DURACAO -------------------------------------
modelo_srpe_lmm <- lmer(
  srpe ~ position * model * day + (1 | athlete) + (1 | session),
  data = dados,
  REML = TRUE
)

# 6. ESTRUTURA DA DURACAO -----------------------------------------------------
modelo_duracao <- lm(duration ~ model * day, data = dados_sessao)

# 7. PSE NUMERICA SEM DURACAO -------------------------------------------------
modelo_rpe_lmm <- lmer(
  rpe ~ position * model * day + (1 | athlete) + (1 | session),
  data = dados,
  REML = TRUE
)

# 8. PSE NUMERICA + DURACAO ---------------------------------------------------
media_duracao_sessoes <- mean(dados_sessao$duration)  # 51 min
dados$duration_c <- dados$duration - media_duracao_sessoes

modelo_rpe_duracao_lmm <- lmer(
  rpe ~ position * model * day + duration_c + (1 | athlete) + (1 | session),
  data = dados,
  REML = TRUE
)

# 9. PSE ORDINAL SEM DURACAO --------------------------------------------------
modelo_rpe_clmm <- clmm(
  rpe_ord ~ position * model * day + (1 | athlete) + (1 | session),
  data = dados,
  link = "logit",
  Hess = TRUE,
  nAGQ = 1
)

# 10. PSE ORDINAL + DURACAO ---------------------------------------------------
modelo_rpe_duracao_clmm <- clmm(
  rpe_ord ~ position * model * day + duration_c + (1 | athlete) + (1 | session),
  data = dados,
  link = "logit",
  Hess = TRUE,
  nAGQ = 1
)

# 11. DIAGNOSTICO AUXILIAR DE ODDS PROPORCIONAIS ------------------------------
# Sem efeitos aleatorios: usado apenas como diagnostico auxiliar.
modelo_clm_diag <- clm(rpe_ord ~ position * model * day, data = dados, link = "logit")
nominal_diag <- nominal_test(modelo_clm_diag)

modelo_clm_dur_diag <- clm(
  rpe_ord ~ position * model * day + duration_c,
  data = dados,
  link = "logit"
)
nominal_dur_diag <- nominal_test(modelo_clm_dur_diag)

# Colapso 8/9/10 apenas para diagnostico de categorias raras
dados$rpe_ord_colapsada <- ordered(ifelse(dados$rpe >= 8, 8, dados$rpe))
modelo_clm_colapsado <- clm(
  rpe_ord_colapsada ~ position * model * day,
  data = dados,
  link = "logit"
)
nominal_colapsado <- nominal_test(modelo_clm_colapsado)

modelo_clm_colapsado_dur <- clm(
  rpe_ord_colapsada ~ position * model * day + duration_c,
  data = dados,
  link = "logit"
)
nominal_colapsado_dur <- nominal_test(modelo_clm_colapsado_dur)

# 12. FUNCOES PARA CONTRASTES -------------------------------------------------
extrair_contrastes <- function(modelo, tipo = c("lmm", "clmm"), duracao = FALSE) {
  tipo <- match.arg(tipo)
  at_arg <- if (duracao) list(duration_c = 0) else NULL

  emm_pos <- emmeans(modelo, ~ position | model * day, at = at_arg)
  emm_mod <- emmeans(modelo, ~ model | position * day, at = at_arg)
  emm_day <- emmeans(modelo, ~ day | position * model, at = at_arg)

  list(
    posicao = as.data.frame(summary(contrast(emm_pos, "pairwise", adjust = "none"), infer = c(TRUE, TRUE))),
    modelo  = as.data.frame(summary(contrast(emm_mod, "pairwise", adjust = "none"), infer = c(TRUE, TRUE))),
    dia     = as.data.frame(summary(contrast(emm_day, "pairwise", adjust = "bonferroni"), infer = c(TRUE, TRUE)))
  )
}

padronizar_contrastes <- function(df, familia, estrutura) {
  inf_col <- if ("lower.CL" %in% names(df)) "lower.CL" else "asymp.LCL"
  sup_col <- if ("upper.CL" %in% names(df)) "upper.CL" else "asymp.UCL"

  if (familia == "Posicao") {
    contexto <- paste(df$model, df$day, sep = " | ")
  } else if (familia == "Modelo") {
    contexto <- paste(df$position, df$day, sep = " | ")
  } else if (familia == "Dia") {
    contexto <- paste(df$position, df$model, sep = " | ")
  } else stop("Familia desconhecida.")

  estimativa <- df$estimate
  ic_inf <- df[[inf_col]]
  ic_sup <- df[[sup_col]]

  out <- data.frame(
    familia = familia,
    contexto = contexto,
    contraste = as.character(df$contrast),
    estrutura = estrutura,
    estimativa = estimativa,
    IC95_inf = ic_inf,
    IC95_sup = ic_sup,
    p = df$p.value,
    direcao = ifelse(estimativa > 0, "Positiva", ifelse(estimativa < 0, "Negativa", "Zero")),
    evidencia = ifelse(ic_inf > 0 | ic_sup < 0, "Sim", "Nao"),
    stringsAsFactors = FALSE
  )
  out$id_contraste <- paste(out$familia, out$contexto, out$contraste, sep = " || ")
  out
}

montar_tabela_estrutura <- function(modelo, estrutura, tipo = "lmm", duracao = FALSE) {
  cts <- extrair_contrastes(modelo, tipo = tipo, duracao = duracao)
  rbind(
    padronizar_contrastes(cts$posicao, "Posicao", estrutura),
    padronizar_contrastes(cts$modelo,  "Modelo",  estrutura),
    padronizar_contrastes(cts$dia,     "Dia",     estrutura)
  )
}

# 13. TABELA-MAE: 24 CONTRASTES x 5 ESTRUTURAS -------------------------------
tabela_srpe <- montar_tabela_estrutura(modelo_srpe_lmm, "PSE_x_duracao")
tabela_rpe <- montar_tabela_estrutura(modelo_rpe_lmm, "PSE_numerica")
tabela_rpe_dur <- montar_tabela_estrutura(modelo_rpe_duracao_lmm, "PSE_numerica_duracao", duracao = TRUE)
tabela_ord <- montar_tabela_estrutura(modelo_rpe_clmm, "PSE_ordinal", tipo = "clmm")
tabela_ord_dur <- montar_tabela_estrutura(modelo_rpe_duracao_clmm, "PSE_ordinal_duracao", tipo = "clmm", duracao = TRUE)

tabela_mae <- rbind(tabela_srpe, tabela_rpe, tabela_rpe_dur, tabela_ord, tabela_ord_dur)
rownames(tabela_mae) <- NULL

# 14. RESUMO DE SENSIBILIDADE ENTRE ESTRUTURAS -------------------------------
ids <- unique(tabela_mae$id_contraste)
resumo_sensibilidade <- do.call(rbind, lapply(ids, function(id) {
  x <- tabela_mae[tabela_mae$id_contraste == id, ]
  data.frame(
    id_contraste = id,
    familia = x$familia[1],
    contexto = x$contexto[1],
    contraste = x$contraste[1],
    direcao_robusta = length(unique(x$direcao)) == 1,
    evidencia_robusta = length(unique(x$evidencia)) == 1,
    n_evidencia = sum(x$evidencia == "Sim"),
    estruturas_com_evidencia = paste(x$estrutura[x$evidencia == "Sim"], collapse = "; "),
    estruturas_sem_evidencia = paste(x$estrutura[x$evidencia == "Nao"], collapse = "; "),
    stringsAsFactors = FALSE
  )
}))
rownames(resumo_sensibilidade) <- NULL

# 15. DECOMPOSICAO DAS DECISOES ANALITICAS -----------------------------------
comparar_estruturas <- function(estrutura_a, estrutura_b) {
  a <- tabela_mae[tabela_mae$estrutura == estrutura_a, c("id_contraste", "direcao", "evidencia")]
  b <- tabela_mae[tabela_mae$estrutura == estrutura_b, c("id_contraste", "direcao", "evidencia")]
  z <- merge(a, b, by = "id_contraste", suffixes = c("_a", "_b"))
  c(
    direcao = sum(z$direcao_a != z$direcao_b),
    evidencia = sum(z$evidencia_a != z$evidencia_b)
  )
}

decomposicao_decisoes <- rbind(
  PSE_x_duracao_vs_PSE = comparar_estruturas("PSE_x_duracao", "PSE_numerica"),
  duracao_aditiva_numerica = comparar_estruturas("PSE_numerica", "PSE_numerica_duracao"),
  PSE_numerica_vs_ordinal = comparar_estruturas("PSE_numerica", "PSE_ordinal"),
  duracao_aditiva_ordinal = comparar_estruturas("PSE_ordinal", "PSE_ordinal_duracao")
)

# 16. DECOMPOSICAO DESCRITIVA PSE x DURACAO ----------------------------------
media_rpe <- aggregate(rpe ~ position + model + day, dados, mean)
media_duracao <- aggregate(duration ~ position + model + day, dados, mean)
media_srpe <- aggregate(srpe ~ position + model + day, dados, mean)

decomposicao <- merge(media_rpe, media_duracao, by = c("position", "model", "day"))
decomposicao <- merge(decomposicao, media_srpe, by = c("position", "model", "day"))
decomposicao$produto_medias <- decomposicao$rpe * decomposicao$duration
decomposicao$componente_covariancia <- decomposicao$srpe - decomposicao$produto_medias
decomposicao$covariancia_percentual <- 100 * decomposicao$componente_covariancia / decomposicao$srpe
decomposicao <- decomposicao[order(decomposicao$model, decomposicao$position, decomposicao$day), ]
rownames(decomposicao) <- NULL

decompor_produto <- function(R_A, t_A, R_B, t_B) {
  data.frame(
    componente_rpe = (R_A - R_B) * t_B,
    componente_duracao = R_B * (t_A - t_B),
    componente_interacao = (R_A - R_B) * (t_A - t_B),
    diferenca_produto = R_A * t_A - R_B * t_B
  )
}

decomp_modelo <- do.call(rbind, lapply(c("Forwards", "Backs"), function(pos) {
  do.call(rbind, lapply(c("StD", "SpD"), function(dia) {
    A <- decomposicao[decomposicao$position == pos & decomposicao$model == "Model_I" & decomposicao$day == dia, ]
    B <- decomposicao[decomposicao$position == pos & decomposicao$model == "Model_II" & decomposicao$day == dia, ]
    cbind(familia = "Modelo", contexto = paste(pos, dia, sep = " | "),
          contraste = "Model_I - Model_II",
          decompor_produto(A$rpe, A$duration, B$rpe, B$duration))
  }))
}))

decomp_dia <- do.call(rbind, lapply(c("Forwards", "Backs"), function(pos) {
  A <- decomposicao[decomposicao$position == pos & decomposicao$model == "Model_II" & decomposicao$day == "StD", ]
  B <- decomposicao[decomposicao$position == pos & decomposicao$model == "Model_II" & decomposicao$day == "EnD", ]
  cbind(familia = "Dia", contexto = paste(pos, "Model_II", sep = " | "),
        contraste = "StD - EnD",
        decompor_produto(A$rpe, A$duration, B$rpe, B$duration))
}))

decomposicao_inversoes <- rbind(decomp_modelo, decomp_dia)
rownames(decomposicao_inversoes) <- NULL

# 17. SENSIBILIDADE: EFEITO ALEATORIO DE SEMANA -------------------------------
modelo_srpe_semana <- lmer(
  srpe ~ position * model * day + (1 | athlete) + (1 | week) + (1 | session),
  data = dados, REML = TRUE
)

modelo_rpe_semana <- lmer(
  rpe ~ position * model * day + (1 | athlete) + (1 | week) + (1 | session),
  data = dados, REML = TRUE
)

tabela_rpe_semana <- montar_tabela_estrutura(modelo_rpe_semana, "PSE_numerica_semana")
principal_rpe <- tabela_mae[tabela_mae$estrutura == "PSE_numerica", ]
comparacao_semana <- merge(
  principal_rpe[, c("id_contraste", "familia", "contexto", "contraste", "estimativa", "direcao", "evidencia")],
  tabela_rpe_semana[, c("id_contraste", "estimativa", "direcao", "evidencia")],
  by = "id_contraste", suffixes = c("_principal", "_semana"), sort = FALSE
)
comparacao_semana$mudou_direcao <- comparacao_semana$direcao_principal != comparacao_semana$direcao_semana
comparacao_semana$mudou_evidencia <- comparacao_semana$evidencia_principal != comparacao_semana$evidencia_semana

# 18. SENSIBILIDADE: OBSERVACOES EXTREMAS DA sRPE -----------------------------
dados$residuo_pad_srpe <- resid(modelo_srpe_lmm) / sigma(modelo_srpe_lmm)
extremos_srpe <- dados[abs(dados$residuo_pad_srpe) > 3, ]
dados_srpe_sem_extremos <- dados[abs(dados$residuo_pad_srpe) <= 3, ]

modelo_srpe_sem_extremos <- lmer(
  srpe ~ position * model * day + (1 | athlete) + (1 | session),
  data = dados_srpe_sem_extremos, REML = TRUE
)

tabela_srpe_sem_extremos <- montar_tabela_estrutura(modelo_srpe_sem_extremos, "PSE_x_duracao_sem_extremos")
principal_srpe <- tabela_mae[tabela_mae$estrutura == "PSE_x_duracao", ]
comparacao_extremos <- merge(
  principal_srpe[, c("id_contraste", "familia", "contexto", "contraste", "estimativa", "direcao", "evidencia")],
  tabela_srpe_sem_extremos[, c("id_contraste", "estimativa", "direcao", "evidencia")],
  by = "id_contraste", suffixes = c("_principal", "_sem_extremos"), sort = FALSE
)
comparacao_extremos$mudou_direcao <- comparacao_extremos$direcao_principal != comparacao_extremos$direcao_sem_extremos
comparacao_extremos$mudou_evidencia <- comparacao_extremos$evidencia_principal != comparacao_extremos$evidencia_sem_extremos

# 19. AUDITORIA AUTOMATICA -----------------------------------------------------
stopifnot(nrow(tabela_mae) == 120)
stopifnot(length(unique(tabela_mae$estrutura)) == 5)
stopifnot(all(table(tabela_mae$estrutura) == 24))
stopifnot(length(unique(tabela_mae$id_contraste)) == 24)
stopifnot(all(table(tabela_mae$id_contraste) == 5))
stopifnot(!anyNA(tabela_mae[c("familia", "contexto", "contraste", "estrutura",
                              "estimativa", "IC95_inf", "IC95_sup", "p", "direcao", "evidencia")]))

evidencia_pelo_ic <- ifelse(tabela_mae$IC95_inf > 0 | tabela_mae$IC95_sup < 0, "Sim", "Nao")
stopifnot(all(evidencia_pelo_ic == tabela_mae$evidencia))

# Resultados congelados esperados da auditoria
stopifnot(nrow(resumo_sensibilidade) == 24)
stopifnot(sum(!resumo_sensibilidade$direcao_robusta & resumo_sensibilidade$evidencia_robusta) == 6)
stopifnot(sum(resumo_sensibilidade$direcao_robusta & !resumo_sensibilidade$evidencia_robusta) == 2)
stopifnot(sum(resumo_sensibilidade$direcao_robusta & resumo_sensibilidade$evidencia_robusta) == 16)
stopifnot(sum(!resumo_sensibilidade$direcao_robusta & !resumo_sensibilidade$evidencia_robusta) == 0)
stopifnot(sum(comparacao_semana$mudou_direcao) == 0)
stopifnot(sum(comparacao_semana$mudou_evidencia) == 1)
stopifnot(nrow(extremos_srpe) == 4)
stopifnot(sum(comparacao_extremos$mudou_direcao) == 1)
stopifnot(sum(comparacao_extremos$mudou_evidencia) == 0)

# 20. EXPORTACAO ---------------------------------------------------------------
write.csv(tabela_mae, "resultados/tabela_mae_contrastes.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(resumo_sensibilidade, "resultados/resumo_sensibilidade_estruturas.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(decomposicao, "resultados/decomposicao_pse_duracao.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(decomposicao_inversoes, "resultados/decomposicao_inversoes.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(comparacao_semana, "resultados/sensibilidade_semana.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(comparacao_extremos, "resultados/sensibilidade_extremos_srpe.csv", row.names = FALSE, fileEncoding = "UTF-8")
write.csv(as.data.frame(decomposicao_decisoes), "resultados/decomposicao_decisoes_analiticas.csv", row.names = TRUE, fileEncoding = "UTF-8")

capture.output(sessionInfo(), file = "resultados/sessionInfo.txt")

cat("\nAnalise concluida com sucesso.\n")
cat("Observacoes:", nrow(dados), "\n")
cat("Contrastes na tabela-mae:", nrow(tabela_mae), "\n")
cat("Contrastes sensiveis a direcao entre estruturas:", sum(!resumo_sensibilidade$direcao_robusta), "\n")
cat("Contrastes sensiveis a evidencia entre estruturas:", sum(!resumo_sensibilidade$evidencia_robusta), "\n")
cat("Mudancas por inclusao de semana (direcao/evidencia):",
    sum(comparacao_semana$mudou_direcao), "/", sum(comparacao_semana$mudou_evidencia), "\n")
cat("Mudancas por exclusao de extremos (direcao/evidencia):",
    sum(comparacao_extremos$mudou_direcao), "/", sum(comparacao_extremos$mudou_evidencia), "\n")
cat("Arquivos salvos em: resultados/\n")