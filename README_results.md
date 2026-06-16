---
title: "Pipeline ITS2 – Comunidades fúngicas de solos da Caatinga (resultados preliminares)"
format: gfm
author: "Daniele Mendes da Silva et al."
---

## 1. Escopo deste README

Este documento resume os **resultados preliminares** obtidos no pipeline ITS2 descrito em `README_methods.qmd`, desde o controle de qualidade dos FASTQ brutos até a filtragem de ASVs fúngicas classificadas e a visualização taxonômica inicial.[file:151][web:66][web:67]  
Aqui são registradas métricas numéricas, flags de qualidade e observações interpretativas que complementam as descrições de métodos.

---

## 2. Organização dos dados brutos e metadados

Os arquivos FASTQ pareados foram organizados em `data/raw/` com subpastas por supersítio (CND, GBR, IBD, JGR), usando o padrão `SAMPLEID_1.fastq` e `SAMPLEID_2.fastq` para leituras forward e reverse, por exemplo `data/raw/CND/CND_06_1.fastq` e `data/raw/CND/CND_06_2.fastq`.[file:135]  
Os arquivos `data/processed/sample-metadata.tsv` e `data/external/metadata.tsv` contêm uma linha por amostra, com identificadores idênticos aos prefixos dos arquivos FASTQ e colunas `Supersitio` e `Subarea` definindo a estrutura espacial (Canudos, Gruta Brejões, Ibiraba_1/IBD1, Ibiraba_2/IBD2, Jaguara).

O manifesto `data/raw/manifest.tsv` foi construído no formato `PairedEndFastqManifestPhred33V2`, com três colunas (`sample-id`, `forward-absolute-filepath`, `reverse-absolute-filepath`), apontando para os caminhos absolutos dos arquivos FASTQ forward e reverse em cada supersítio.[file:135][web:79]  
Esse manifest é usado diretamente pelo QIIME2 para importar as leituras demultiplexadas, garantindo consistência entre os IDs de amostra dos metadados e dos arquivos de sequência.

---

## 3. Controle de qualidade dos FASTQ (FastQC + MultiQC)

### 3.1. Configuração e ferramentas

O controle de qualidade inicial foi realizado com FastQC (v0.12.1), executado sobre todos os arquivos FASTQ em `data/raw/*/*.fastq`, e os relatórios individuais foram agregados em um relatório MultiQC único (`results/multiqc/multiqc_report.html`).[file:151][web:66][web:67]  
Essa etapa foi executada dentro do ambiente Conda `qiime2-amplicon-2024.10`, no qual FastQC e MultiQC foram instalados adicionalmente ao QIIME2 conforme descrito em `README_methods.qmd`.[web:136][web:145]

O relatório MultiQC inclui, para cada amostra e direção de leitura, os módulos padrão do FastQC: per base sequence quality, per sequence GC content, sequence length distribution, sequence duplication levels, overrepresented sequences, adapter content e sumário de status (heatmap de PASS/WARN/FAIL).[file:151][web:66]  

### 3.2. Sequences duplicadas e overrepresented sequences

O relatório MultiQC mostra que as bibliotecas apresentam **altos níveis de sequências duplicadas e de sequências sobre‑representadas**, com múltiplas sequências específicas aparecendo em todas ou quase todas as amostras, como esperado em dados de amplicon ITS2 com priming altamente conservado.[file:151][file:7]  
A tabela “Top overrepresented sequences” indica que a sequência mais abundante (por exemplo, `GAAATGCGATAAGTAATGTGAATTGCAGAATTCAGTGAATCATCGAATCT`) aparece em todas as 50 amostras, representando ≈10,5% de todas as leituras, seguida por outras variantes próximas com frequências entre ~10% e <1% do total de reads.[file:151]

Outras sequências ITS‑like altamente frequentes também aparecem em dezenas de amostras, refletindo a natureza direcionada do metabarcoding (amplicons de alta cópia, primers específicos) e não um problema de artefato de sequenciamento aleatório como seria interpretado em bibliotecas de shotgun.[file:151][file:7]  
Esses níveis elevados de duplicação e de sequências sobre‑representadas são esperados para bibliotecas de amplicon e serão tratados explicitamente pelo DADA2 na etapa de inferência de ASVs, que modela a distribuição de erros e colapsa reads idênticas em variantes únicas com suporte de abundância.[file:7][web:49]

### 3.3. Conteúdo de adaptadores

Na seção “Adapter Content”, o MultiQC sumariza a presença acumulativa de adaptadores ao longo do comprimento das leituras, reportando apenas amostras com ≥0,1% de contaminação por adaptadores.[file:151]  
Os perfis indicam que há detecção de sequências associadas aos adaptadores de Illumina e/ou primers, principalmente em regiões mais ao final dos reads, o que reforça a necessidade da etapa subsequente de remoção explícita dos primers ITS3/ITS4 com o plugin Cutadapt do QIIME2.[file:151][web:81]

O padrão de aumento cumulativo de adaptadores ao longo dos ciclos é compatível com a descrição do próprio FastQC, em que uma vez que a sequência de adaptador é detectada numa leitura, ela é contada até o final do read, levando a curvas ascendentes nas posições finais.[file:151][web:66]  
Esse comportamento será corrigido no pipeline ao aplicar o comando `qiime cutadapt trim-paired` com as sequências exatas dos primers ITS3 e ITS4, descartando leituras onde os primers não forem encontrados (`--p-discard-untrimmed`), conforme descrito no README de métodos.[web:81][file:7]

### 3.4. Interpretação geral do QC

O heatmap de “Status Checks” do MultiQC resume os módulos do FastQC em avaliações do tipo PASS/WARN/FAIL para cada amostra, mostrando o padrão típico de bibliotecas de amplicon, com alertas concentrados em módulos de **sequence duplication levels**, **overrepresented sequences** e, em alguns casos, **adapter content**.[file:151][web:66]  
Tais alerts não indicam necessariamente problemas de qualidade para enriquecimento por amplicon, mas sim a alta redundância de reads e a forte presença de sequências alvo, que serão tratadas pelas etapas de removal de primers (Cutadapt) e inferência de ASVs (DADA2).[file:7][web:49]

Globalmente, o QC confirma que as bibliotecas de ITS2 de solo da Caatinga são adequadas para prosseguimento no pipeline, desde que se mantenham as etapas planejadas de remoção de primers/adaptadores e de denoising, e que a interpretação dos módulos de FastQC seja feita no contexto de metabarcoding em vez de bibliotecas de shotgun randomizadas.[file:151][file:7]  

---

## 4. Processamento no QIIME2

### Profundidade de sequenciamento e controle de qualidade inicial

As leituras pareadas ITS foram importadas no QIIME 2 (versão 2024.10) por meio de um arquivo de manifest no formato PairedEndFastqManifestPhred33V2, gerando o artefato `demux.qza` com as sequências demultiplexadas.[web:152][web:194] A etapa `qiime demux summarize` produziu o sumário `demux.qzv`, a partir do qual foi avaliada a distribuição de profundidade de sequenciamento por amostra e os perfis de qualidade das leituras forward e reverse.[web:194][web:197]

Foram obtidas 50 amostras com leituras pareadas, totalizando 9.786.539 reads forward e 9.786.539 reads reverse antes de qualquer filtragem de qualidade.[web:194] O número de reads por amostra variou de 87.112 a 219.097, com mediana de 204.585 reads por amostra para ambas as direções, indicando profundidade de sequenciamento elevada e relativamente uniforme entre as amostras.[web:194] Esses valores atendem com folga as recomendações de profundidade para estudos de metabarcoding fúngico, garantindo poder estatístico adequado para as análises de diversidade e composição de comunidades.[file:170][file:172][web:200]

A inspeção das curvas de qualidade no sumário `demux.qzv` mostrou perfis de qualidade bastante homogêneos entre amostras para ambas as direções.[web:194][web:197] As distribuições de comprimento foram estreitas, com medianas de 225 nucleotídeos para as leituras forward e 226 nucleotídeos para as reverse (intervalos interquartis aproximadamente 223–228 nucleotídeos em ambas as direções).[web:194] A mediana dos escores de qualidade Phred permaneceu elevada (≈ Q37–38) até cerca de 200–210 pares de bases, com aumento gradual da dispersão apenas nos ciclos finais, mais pronunciado nas leituras reverse, padrão compatível com dados MiSeq de amplicons ITS.[web:194][web:195][web:204] Esses perfis orientaram a escolha de comprimentos de truncamento em torno de 220 pb (forward) e 210 pb (reverse) na etapa de denoising com DADA2, de forma a remover as regiões de menor qualidade preservando sobreposição suficiente entre as leituras pareadas.[web:195][web:204][web:200]

### 4.1. Desempenho do DADA2

Após o denoising com DADA2, a maior parte das leituras manteve-se ao longo do pipeline, com aproximadamente 72–94% das reads por amostra passando na etapa de filtragem inicial e 80–91% dos pares de leituras sendo corretamente unidos.[web:195][web:204] A proporção final de reads não quiméricas variou de cerca de 50% a 88% do total de reads de entrada por amostra, com a maioria das amostras retendo em torno de 70–80% das leituras originais, desempenho comparável ao observado em estudos recentes de metabarcoding fúngico em solos.[file:170][file:171][web:200]

### 4.2. Filtragem taxonômica

A filtragem taxonômica e de abundância resultou em uma tabela final contendo 6.664 ASVs atribuídas ao Reino Fungi e 6.432.099 leituras distribuídas em 50 amostras.[web:194][web:219] A profundidade de sequenciamento por amostra após o filtro variou de 55.394 a 176.426 reads, com mediana de 137.208 reads por amostra, indicando que a remoção de ASVs não fúngicas, não classificadas e de baixa frequência (frequência total < 5 reads) manteve elevada a cobertura de sequenciamento por amostra.[web:194][file:169][file:172] A distribuição das frequências por feature apresentou mediana de 57 reads e máximo de 332.916 reads por ASV, sugerindo que a maioria das variantes retidas possui suporte de abundância suficiente para análises subsequentes de diversidade e composição de comunidades.[web:219][file:171]

### 4.3. Visualização taxonômica

A visualização interativa `results/qiime2/taxa-barplot_fungi.qzv` foi gerada a partir da tabela filtrada (`workflow/04-filtering/table_fungi_classified.qza`), da taxonomia UNITE (`workflow/03-taxonomy/taxonomy.qza`) e dos metadados em `data/external/metadata.tsv`. Esse arquivo permite explorar a composição relativa de táxons fúngicos por amostra, supersítio e subárea.

## 5. Estado atual e próximas etapas

O conjunto de dados processado está pronto para análises ecológicas downstream em R ou QIIME2, usando como entrada principal `workflow/04-filtering/table_fungi_classified.qza`, `workflow/04-filtering/rep-seqs_fungi_classified.qza`, `workflow/03-taxonomy/taxonomy.qza` e `data/external/metadata.tsv`.

As próximas etapas recomendadas são:

- Exportar as tabelas finais para `results/tables/` em formatos tabulares (`.tsv`/`.biom`) para análise em R.
- Calcular diversidade alfa e beta por supersítio/subárea, com escolha explícita de profundidade de rarefação ou alternativa sem rarefação.
- Gerar ordenações, testes estatísticos e figuras finais em `results/figures/`.
- Versionar os scripts ainda ausentes para DADA2 e classificação taxonômica, garantindo rastreabilidade completa entre comandos documentados e arquivos em `scripts/`.

### Figura 3 – desempenho da etapa DADA2

A Figura 3 resume o desempenho da etapa de denoising com DADA2 para os amplicons ITS2 de fungos de solos da Caatinga. Os dados foram obtidos a partir do artefato `workflow/02-dada2/denoising-stats.qza`, exportado com `qiime tools export` para `data/processed/qc_export/denoising_stats_export/stats.tsv`, e combinados com os metadados amostrais (`data/external/metadata.tsv`) em R.[web:27][web:155][web:182]

A Figura 3A mostra a distribuição do número de reads por amostra antes da inferência de ASVs (coluna `input`) e após a remoção de quimeras (`non-chimeric`). A Figura 3B apresenta a porcentagem de reads retidas em relação ao total de entrada em três momentos do pipeline do DADA2: após a filtragem de qualidade (`percentage of input passed filter`), após a junção de pares (`percentage of input merged`) e após a remoção de quimeras (`percentage of input non-chimeric`).[web:155][web:153]

Os gráficos foram gerados em R com `ggplot2` a partir do script `scripts/fig3_qc_dada2.R`, e exportados para `results/figures/Fig3A_QC_DADA2_reads.(png|tiff)` e `results/figures/Fig3B_QC_DADA2_retention.(png|tiff)` em 600 dpi, formato adequado para submissão em periódicos.[web:186][web:189]


## Legenda
Figure 3A. Number of reads per sample before DADA2 denoising and after removal of chimeric sequences. Values were obtained from the denoising-stats.qza artifact exported as stats.tsv

## Composição taxonômica das comunidades fúngicas

A composição taxonômica das comunidades fúngicas mostrou clara predominância do filo **Ascomycota** em todos os supersítios avaliados, seguido por **Basidiomycota**, enquanto os demais filos apresentaram baixa abundância relativa [file:247]. Entre os sítios, Ibiraba apresentou a maior contribuição relativa de Basidiomycota, mas o padrão geral permaneceu fortemente dominado por Ascomycota em todas as localidades [file:247].

Em níveis taxonômicos mais refinados, observou-se predominância da classe **Eurotiomycetes** em todos os supersítios, com variações locais na participação de **Sordariomycetes**, **Dothideomycetes** e **Agaricomycetes** [file:248]. Gruta_Brejoes destacou-se por maior abundância relativa de Sordariomycetes, enquanto Ibiraba e Jaguara apresentaram maior domínio de Eurotiomycetes [file:248].

Esses resultados indicam que, embora os supersítios compartilhem uma estrutura taxonômica geral semelhante, há diferenças na contribuição relativa de grupos específicos, sobretudo em nível de classe, o que sugere heterogeneidade local na composição das comunidades fúngicas [file:247][file:248].

::: {#fig-taxonomic-composition layout-ncol=2}
![Composição taxonômica em nível de filo.](/path/to/Fig4A_taxonomic_composition_phylum.jpg)

![Composição taxonômica em nível de classe.](/path/to/Fig4B_taxonomic_composition_class-2.jpg)

**Figura 4.** Composição taxonômica relativa das comunidades fúngicas nos supersítios Canudos, Gruta_Brejoes, Ibiraba e Jaguara. **(A)** Abundância relativa dos principais filos identificados, evidenciando predominância de Ascomycota em todos os sítios. **(B)** Abundância relativa das principais classes taxonômicas, com destaque para Eurotiomycetes e variações na participação de Sordariomycetes, Dothideomycetes e Agaricomycetes entre os supersítios. Táxons pouco abundantes foram agrupados em “Others”. [file:247][file:248]
:::

Figure 4A. Taxonomic composition at the phylum level of soil fungal communities in the supersites Canudos, Gruta_Brejoes, Ibiraba, and Jaguara, expressed as relative read abundance. The stacked bar plot shows a strong dominance of Ascomycota across all sites, with a secondary contribution of Basidiomycota and low representation of other phyla (e.g. Glomeromycota, Mucoromycota, Rozellomycota, Chytridiomycota) and unclassified groups, which are aggregated into the “Others” category.

Figure 4B. Taxonomic composition at the class level of soil fungal communities in the supersites Canudos, Gruta_Brejoes, Ibiraba, and Jaguara, shown as relative abundance. Eurotiomycetes dominates most samples, followed by variable contributions of Sordariomycetes, Dothideomycetes, Agaricomycetes, and other classes (including Glomeromycetes, Leotiomycetes, Orbiliomycetes, and incertae sedis groups), while rare taxa are grouped as “Others” and unassigned reads as “Unclassified


Riqueza (nº de ASVs) por supersítio

riqueza_supersitio
      Canudos Gruta_Brejoes       Ibiraba       Jaguara 
         2771           947          2553          1391 

exclusivos_por_supersitio
      Canudos Gruta_Brejoes       Ibiraba       Jaguara 
         2271           782          2087          1136 


## Resultados gerados até o momento

Até esta etapa, foi estruturado um conjunto de tabelas e figuras derivadas da matriz de abundância de ASVs fúngicos, da taxonomia expandida e do arquivo de metadata das amostras, permitindo análises em nível de amostra e em nível de supersítio.[cite:738][cite:739]

### 1. Integração entre abundância e taxonomia

A partir da tabela de abundância de fungos (`asv_tab_fungi`) e da tabela taxonômica expandida (`tax_tidy`), foi realizada a interseção dos identificadores de ASVs para manter apenas os registros comuns entre os dois objetos. Esse procedimento garante que cada ASV presente na matriz de abundância possua uma anotação taxonômica correspondente.[cite:738]

Após o alinhamento, foi gerada a tabela consolidada:

```r
results/asv_table_fungi_with_taxonomy.tsv
```

Esse arquivo reúne, em uma única estrutura, os identificadores de ASVs, a taxonomia expandida e a abundância de reads por amostra. Ele constitui a tabela-base para análises ecológicas por amostra, incluindo diversidade alfa, diversidade beta e descrições taxonômicas em diferentes níveis hierárquicos.[cite:738][cite:739]

### 2. Integração com o metadata e agregação por supersítio

O arquivo `metadata.tsv`, localizado na raiz do projeto, foi utilizado para associar cada amostra ao seu respectivo supersítio, seguindo o formato tabular amplamente empregado em fluxos de trabalho do QIIME 2.[cite:563][cite:574]

Com base nessa associação, a abundância das ASVs foi agregada por supersítio, somando-se as contagens de todas as amostras pertencentes ao mesmo grupo. Como resultado, foi gerada a tabela:

```r
results/asv_por_supersitio_com_taxonomia.tsv
```

Essa tabela permite explorar a composição fúngica agregada por supersítio, sendo útil para análises descritivas, construção de gráficos de composição taxonômica e comparação geral entre áreas.[cite:739]

### 3. Tabela binária de presença e ausência

A matriz agregada por supersítio também foi convertida para formato binário, em que valores maiores que zero foram recodificados como presença (1) e valores iguais a zero como ausência (0). Esse procedimento gerou o arquivo:

```r
results/asv_por_supersitio_presenca_ausencia.tsv
```

Essa estrutura é apropriada para análises focadas em ocorrência de ASVs, como contagem de ASVs exclusivas e compartilhadas entre supersítios, além de métricas de dissimilaridade baseadas em presença/ausência.[cite:665][cite:645]

### 4. Diversidade alfa por amostra

A diversidade alfa foi calculada em nível de amostra a partir da tabela `results/asv_table_fungi_with_taxonomy.tsv`, após separação entre as colunas taxonômicas e a matriz de abundância. A matriz foi transposta para o formato amostras x ASVs, conforme requerido pelo pacote `vegan` para cálculos ecológicos.[cite:659][cite:671]

Foram calculadas três métricas por amostra:

- riqueza observada (`specnumber`), correspondente ao número de ASVs observadas por amostra.[cite:659]
- índice de Shannon (`diversity(..., index = "shannon")`), que sintetiza riqueza e equabilidade.[cite:736][cite:741]
- índice de Simpson (`diversity(..., index = "simpson")`), como medida adicional de diversidade intra-amostra.[cite:671]

Os resultados foram organizados na tabela:

```r
results/alpha_diversidade_por_amostra.tsv
```

Esse arquivo contém uma linha por amostra, com os índices de diversidade alfa e as variáveis do metadata associadas, possibilitando análises estatísticas e comparações entre grupos amostrais.[cite:659][cite:737]

### 5. Figuras de diversidade alfa

Com base na tabela `alpha_diversidade_por_amostra.tsv`, foram preparados boxplots para comparar a diversidade alfa entre supersítios, com rótulos em inglês para facilitar o uso em manuscritos e apresentações.[cite:659]

As figuras planejadas/salvas foram:

```r
results/figures/alpha_shannon_by_supersite.png
results/figures/alpha_shannon_by_supersite.pdf
results/figures/alpha_richness_by_supersite.png
results/figures/alpha_richness_by_supersite.pdf
```

A exportação foi feita com `ggplot2::ggsave`, que permite salvar os gráficos em múltiplos formatos com definição explícita de dimensões e resolução.[cite:708]

### 6. Uso analítico dos arquivos produzidos

Os produtos gerados até aqui já permitem diferentes frentes de análise:[cite:739][cite:645]

- `asv_table_fungi_with_taxonomy.tsv`: base para diversidade alfa e beta por amostra, ordenações, testes multivariados e descrições taxonômicas detalhadas.[cite:737][cite:739]
- `asv_por_supersitio_com_taxonomia.tsv`: base para composição por phylum, class ou outros níveis taxonômicos em escala de supersítio.[cite:739]
- `asv_por_supersitio_presenca_ausencia.tsv`: base para análises de exclusividade, compartilhamento e dissimilaridade por ocorrência.[cite:665][cite:645]
- `alpha_diversidade_por_amostra.tsv`: base para boxplots, testes estatísticos entre grupos e integração com variáveis ambientais ou experimentais.[cite:736][cite:741]

### 7. Próximas análises recomendadas

Com a estrutura atual, os próximos passos mais diretos incluem a execução de testes estatísticos para comparar diversidade alfa entre supersítios, a geração de gráficos de composição taxonômica por phylum e class, e a análise de diversidade beta por métricas como Bray-Curtis e Jaccard.[cite:737][cite:665][cite:645]

## Resultados da diversidade beta por Bray-Curtis

A análise de diversidade beta baseada em Bray-Curtis indicou estruturação clara das comunidades fúngicas do solo em função dos supersítios, com separação visual consistente entre grupos no espaço de ordenação e efeito estatisticamente significativo do fator supersítio sobre a composição das amostras.[cite:647][cite:777]

### 1. Qualidade da ordenação NMDS

A ordenação por NMDS apresentou valor de *stress* igual a 0.1164705, indicando ajuste aceitável da representação bidimensional para descrever o padrão multivariado de dissimilaridade entre amostras. Em termos práticos, esse valor sugere que a projeção em dois eixos preserva de forma satisfatória a estrutura das distâncias Bray-Curtis observadas no conjunto de dados.[cite:781][cite:782]

### 2. Padrão visual de separação entre supersítios

No gráfico NMDS, as amostras formaram agrupamentos compatíveis com a origem por supersítio, evidenciando diferenciação composicional entre as comunidades fúngicas avaliadas. O padrão visual mostra que os grupos ocupam regiões distintas do espaço de ordenação, com separação particularmente evidente entre conjuntos como Ibiraba, Gruta_Brejoes, Canudos e Jaguara.[cite:647]

[cite:784]

A distribuição dos pontos sugere que os supersítios compartilham apenas parcialmente sua composição fúngica, reforçando a hipótese de que condições ambientais locais contribuem para a organização da comunidade em escala espacial.[cite:645][cite:665]

### 3. PERMANOVA por supersítio

A PERMANOVA revelou efeito significativo do fator supersítio sobre a composição fúngica das amostras (`F = 7.5542`, `R2 = 0.33006`, `p = 0.001`). Esse resultado indica que as diferenças entre supersítios explicam aproximadamente 33.0% da variação total observada na matriz de dissimilaridade Bray-Curtis.[cite:777]

Os valores obtidos foram:

```r
Df = 3
SumOfSqs = 6.7518
R2 = 0.33006
F = 7.5542
Pr(>F) = 0.001
```

O componente residual correspondeu a `R2 = 0.66994`, mostrando que parte substancial da variação ainda permanece associada a diferenças intra-supersítio ou a outros fatores não incluídos neste modelo inicial.[cite:777][cite:645]

### 4. Interpretação ecológica inicial

Em conjunto, o NMDS e a PERMANOVA indicam que os supersítios não representam apenas agrupamentos geográficos, mas também unidades com composição fúngica distinta. Esse padrão é coerente com a expectativa de que gradientes locais de solo, vegetação, microclima ou histórico ambiental influenciem a montagem das comunidades em ecossistemas secos.[cite:645][cite:665]

Como o efeito de supersítio foi significativo e apresentou magnitude moderada (`R2 ≈ 0.33`), essa variável deve ser mantida como um eixo central nas análises subsequentes, incluindo composição taxonômica, diversidade alfa e comparações multivariadas adicionais.[cite:777][cite:736]

### 5. Próximos passos recomendados

A partir desses resultados, os próximos passos mais indicados são:[cite:647][cite:777]

- testar diversidade beta também com Jaccard para verificar se o padrão se mantém em presença/ausência.[cite:647][cite:761]
- avaliar dispersão multivariada entre grupos para complementar a interpretação da PERMANOVA.[cite:777]
- integrar essas diferenças composicionais com os resultados de diversidade alfa e com a composição por phylum e class já estruturadas no pipeline.[cite:736][cite:659]

## Resultados da dispersão multivariada por Bray-Curtis

A análise de homogeneidade de dispersões multivariadas indicou que os supersítios diferem significativamente em sua variabilidade interna quando avaliados pela matriz de dissimilaridade Bray-Curtis. Esse resultado mostra que, além das diferenças entre centróides observadas na PERMANOVA, há também diferenças na dispersão das amostras dentro dos grupos.[cite:788][cite:807]

### 1. Distância média ao centro do grupo

As distâncias médias à mediana do grupo foram distintas entre os supersítios, com os seguintes valores observados:[cite:788]

```r
Canudos        = 0.4791
Gruta_Brejoes  = 0.5983
Ibiraba        = 0.5212
Jaguara        = 0.3953
```

Entre os grupos avaliados, `Gruta_Brejoes` apresentou a maior dispersão média, enquanto `Jaguara` apresentou a menor. Esse padrão sugere que a comunidade fúngica em Gruta_Brejoes é mais heterogênea entre amostras, ao passo que Jaguara apresenta composição mais internamente consistente.[cite:788][cite:803]

### 2. Teste PERMDISP

O teste permutacional de homogeneidade de dispersões foi significativo (`F = 13.802`, `p = 0.001`), indicando rejeição da hipótese nula de dispersão homogênea entre supersítios.[cite:789][cite:799]

Os valores obtidos foram:

```r
Df = 3
Sum Sq = 0.17097
Mean Sq = 0.056988
F = 13.802
Pr(>F) = 0.001
```

A significância observada demonstra que os supersítios não diferem apenas na posição média de suas comunidades no espaço multivariado, mas também na amplitude de variação interna entre suas amostras.[cite:788][cite:807]

### 3. Implicação para a interpretação da PERMANOVA

Como a PERMANOVA anterior também foi significativa, a presença de PERMDISP significativo indica que a interpretação das diferenças entre supersítios deve ser feita com cautela. Em termos práticos, parte da separação detectada entre grupos pode refletir não apenas diferenças na composição média, mas também diferenças na heterogeneidade interna das comunidades.[cite:801][cite:805]

Isso não invalida o resultado da PERMANOVA, mas indica que o efeito de supersítio sobre a beta diversidade é composto por pelo menos dois componentes: deslocamento entre grupos e diferença de dispersão entre grupos. Por isso, a evidência atual sustenta que os supersítios são distintos, mas a discussão deve reconhecer explicitamente que a variabilidade interna não é homogênea entre eles.[cite:807][cite:789]

### 4. Leitura ecológica inicial

A maior dispersão observada em `Gruta_Brejoes` pode refletir maior heterogeneidade ambiental local, maior variação entre micro-hábitats ou maior sensibilidade da comunidade fúngica a fatores finos de solo e vegetação. Em contraste, a menor dispersão em `Jaguara` sugere uma comunidade relativamente mais uniforme entre as amostras desse supersítio.[cite:803][cite:807]

Esse padrão reforça a necessidade de integrar a interpretação da beta diversidade com resultados de composição taxonômica, diversidade alfa e, idealmente, variáveis ambientais medidas em cada área.[cite:807][cite:647]

### 5. Próximos passos recomendados

Com base nesses resultados, os próximos passos mais úteis são:[cite:789][cite:807]

- executar a análise com Jaccard para verificar se o padrão de heterogeneidade também aparece em presença/ausência.[cite:759]
- realizar comparações par-a-par entre supersítios, se necessário, usando recursos adicionais do `betadisper` para identificar quais grupos diferem mais fortemente em dispersão.[cite:789][cite:802]
- discutir conjuntamente NMDS, PERMANOVA e PERMDISP no texto de resultados, deixando explícito que houve diferença entre grupos e também heterogeneidade interna desigual.[cite:801][cite:805]

## NMDS, PERMANOVA e betadisper por supersítio

A beta-diversidade fúngica foi avaliada com base em uma matriz binária de presença/ausência derivada da tabela de abundância de ASVs. Essa transformação permitiu que a dissimilaridade entre amostras fosse estimada exclusivamente em termos de ocorrência e substituição de ASVs, sem influência da magnitude de abundância. [web:911]

A ordenação NMDS baseada na distância de Jaccard binária apresentou baixo valor de *stress* (stress = 0,0815), indicando bom ajuste entre a configuração bidimensional e a matriz original de dissimilaridades. Isso sugere que a posição relativa das amostras no plano de ordenação representa adequadamente a estrutura composicional da comunidade fúngica. [web:772][web:911]

As 50 amostras analisadas distribuíram-se entre quatro supersítios: Canudos (n = 12), Gruta_Brejoes (n = 12), Ibiraba (n = 21) e Jaguara (n = 5). A distribuição das amostras no espaço de ordenação evidenciou estruturação da comunidade de acordo com a identidade do supersítio, com agrupamentos visíveis na ordenação NMDS. [web:911]

A PERMANOVA confirmou que a composição da comunidade fúngica diferiu significativamente entre supersítios (F = 4,283; R² = 0,218; p = 0,001). O fator `Supersitio` explicou aproximadamente 21,8% da variação total observada na matriz de dissimilaridade de Jaccard binária, indicando que a identidade geográfica dos supersítios constitui um componente importante da estrutura da comunidade fúngica. [web:913][web:914]

A análise de homogeneidade das dispersões multivariadas com `betadisper` indicou diferenças significativas na dispersão entre supersítios (ANOVA sobre as distâncias ao centróide: F = 31,364; p = 0,001). O teste de permutação associado (`permutest`) confirmou esse resultado (F = 31,364; p = 0,001; 999 permutações), sugerindo que parte da separação observada entre supersítios também reflete diferenças na heterogeneidade interna das comunidades fúngicas em cada grupo. [web:877][web:800][web:916]

As figuras correspondentes foram salvas em `results/nmds_permanova/` nos formatos PNG e PDF, incluindo o gráfico de NMDS (`nmds_jaccard_supersitio.png`/`.pdf`) e o gráfico de dispersão multivariada (`betadisper_supersitio.png`/`.pdf`). Os arquivos auxiliares exportados também incluem a matriz binária por amostra, a matriz de distâncias de Jaccard, os escores do NMDS, a tabela da PERMANOVA e as tabelas do `betadisper`, garantindo a reprodutibilidade completa da análise. [web:910][web:915]