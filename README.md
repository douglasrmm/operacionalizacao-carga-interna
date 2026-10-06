# A operacionalização matemática pode afetar inferências sobre a Carga Interna de Treinamento

Repositório de materiais para reprodução das análises do Trabalho de Conclusão de Curso (TCC) de MBA em Data Science & Analytics.

O estudo investigou a sensibilidade das inferências sobre a Carga Interna de Treinamento (CIT) a diferentes formas de representar a Percepção Subjetiva de Esforço da sessão (PSEs) e sua relação com a duração da sessão.

## Conteúdo

- `analise_tcc.R`: script final das análises estatísticas.
- `hu_rpe_analitico.csv`: banco analítico utilizado pelo script.
- `resultados/`: pasta criada automaticamente pelo script, contendo as tabelas de resultados e o `sessionInfo.txt`.

## Reprodução das análises

O script foi estruturado para ser executado com `analise_tcc.R` e `hu_rpe_analitico.csv` no mesmo diretório de trabalho.

Pacotes necessários no R:

```r
install.packages(c("lme4", "lmerTest", "ordinal", "emmeans", "pbkrtest"))
```

Em seguida, execute:

```r
source("analise_tcc.R")
```

O script realiza auditorias automáticas do banco e dos principais resultados e cria a pasta `resultados/` com os arquivos exportados.

## Dados e atribuição

Os dados utilizados são provenientes de:

Hu, X.; Boisbluche, S.; Philippe, K.; Maurelli, O.; Ren, X.; Li, S.; Xu, B.; Prioux, J. (2024). *Position-specific workload of professional rugby union players during tactical periodization training*. **PLOS ONE, 19**(3), e0288345. https://doi.org/10.1371/journal.pone.0288345

O artigo e seus dados foram disponibilizados sob a licença **Creative Commons Attribution 4.0 International (CC BY 4.0)**, que permite uso, distribuição e reprodução mediante atribuição apropriada aos autores e à fonte original.

O arquivo `hu_rpe_analitico.csv` corresponde à versão analítica preparada para as análises deste TCC. A autoria e a fonte dos dados originais permanecem atribuídas a Hu et al. (2024).

## Estruturas analíticas

Foram comparadas cinco estruturas:

1. PSEs × duração;
2. PSEs numérica;
3. PSEs numérica com duração como covariável;
4. PSEs ordinal;
5. PSEs ordinal com duração como covariável.

As análises incluíram modelos lineares mistos, modelos de ligação cumulativa mistos, contrastes pós-estimação, decomposição da operacionalização multiplicativa e análises de sensibilidade.

## Autor

Douglas Miranda
