---
title: "Pipeline ITS2 – Comunidades fúngicas de solos da Caatinga (métodos)"
format: gfm
author: "Daniele Mendes da Silva et al."
---

## 1. Visão geral

Este documento descreve, de forma reprodutível e alinhada aos princípios FAIR, o pipeline de análise da região ITS2 do rDNA de fungos de solos de quatro supersítios de Caatinga na Bahia (Canudos, Gruta Brejões, Ibiraba e Jaguara).  
O fluxo vai dos arquivos FASTQ brutos até a geração de artefatos QIIME2 contendo ASVs ITS2 anotadas taxonomicamente, prontos para análises ecológicas downstream em R.[file:7][file:6]

Os dados de sequência são amplicons da região ITS2 gerados em plataforma Illumina MiSeq (2×300 ciclos, kit v3), com primers ITS3/ITS4, em serviço de sequenciamento terceirizado.[file:5][file:6]  
Os comandos principais são mantidos em scripts versionados na pasta `scripts/` e os artefatos intermediários em `workflow/` e `results/`, garantindo rastreabilidade desde os FASTQ até os artefatos `.qza`/`.qzv`.[cite:122][file:4]

---

## 2. Dados de entrada e metadados

### 2.1. Estrutura de diretórios

A raiz do projeto é:

`/home/user/Documentos/1_Projetos/2026_UEFS_Projeto_metabarcoding_fungos_2026/02_Work/fungos_caatinga_V2`

A estrutura lógica de pastas utilizada é:

```bash
data/
  raw/                     # FASTQ brutos e manifest
  external/                # metadados usados nos comandos QIIME2 e bancos UNITE
  processed/               # manifest original, metadados processados e planilha FAIR
env/
  qiime2-amplicon-2024.10-py310-linux-conda.yml
workflow/
  01-import/
  01b-cutadapt/
  02-dada2/
  03-taxonomy/
  04-filtering/
scripts/
  00-fastqc_multiqc.sh
  01-import.sh
  01b-cutadapt.sh
  04-filtering.sh
  05-taxa-barplot.sh
results/
  fastqc/
  multiqc/
  qiime2/
  tables/
  figures/
README_methods.qmd
README_results.qmd
```

Observação de reprodutibilidade: os artefatos de DADA2 (`workflow/02-dada2/`) e taxonomia (`workflow/03-taxonomy/`) estão presentes, mas os scripts correspondentes (`02-dada2.sh` e `03-taxonomy.sh`) ainda não estão versionados em `scripts/`. Os comandos usados nessas etapas estão documentados nas seções 4.3 e 4.4.

Os arquivos FASTQ estão organizados em subpastas por supersítio dentro de `data/raw/` (CND, GBR, IBD, JGR), com nomes no formato `SAMPLEID_1.fastq` (forward) e `SAMPLEID_2.fastq` (reverse), por exemplo:

```text
data/raw/CND/CND_01_1.fastq
data/raw/CND/CND_01_2.fastq
...
```

### 2.2. Metadados QIIME2 (`data/processed/sample-metadata.tsv` e `data/external/metadata.tsv`)

Os arquivos `data/processed/sample-metadata.tsv` e `data/external/metadata.tsv` registram os metadados usados no pipeline. Ambos contêm uma linha por amostra, com identificadores consistentes com os FASTQ e com o manifest. O arquivo usado nos comandos finais de visualização taxonômica é `data/external/metadata.tsv`, com as colunas:

- `SampleID` – identificador único de cada amostra (ex.: CND_01, CND_03, ..., GBR_01, IBD1_05, IBD2_14, JGR_06).  
- `Supersitio` – Canudos, Gruta_Brejoes, Ibiraba, Jaguara.  
- `Subarea` – para Ibiraba: Ibiraba_1 (IBD1) ou Ibiraba_2 (IBD2); para os demais, igual ao supersítio.[cite:126][cite:129]

Os valores de `SampleID`/`sample-id` são consistentes com os prefixos dos arquivos FASTQ em `data/raw/` (para cada amostra existe um par `sample-id_1.fastq` / `sample-id_2.fastq`).[cite:126]

**Campos adicionais planejados (a serem adicionados depois):**

- `collection_date` – data de coleta (YYYY-MM-DD).  
- `depth` – profundidade (cm).  
- `lat_lon` – coordenadas geográficas.  
- `vegetation_type` / `land_use` – tipo de vegetação/uso do solo.  
- Atributos físico‑químicos de solo (pH, textura, matéria orgânica etc.).[file:6][file:8]

### 2.3. Metadados MIxS/MIMARKS (`data/processed/metadata_fair_ITS2.xlsx`)

O arquivo `data/processed/metadata_fair_ITS2.xlsx` concentra a planilha FAIR/MIxS-MIMARKS para amostras de solo (`MIMARKS.survey` + pacote Soil), com uma linha por amostra (`samp_name` = `sample-id`).[web:24][web:25]

Campos principais incluídos:

- `project_name = Fungos_solo_Caatinga`  
- `env_broad_scale`, `env_local_scale`, `env_medium` – termos ENVO para Caatinga e solo (ex.: bioma Caatinga, solo arenoso/latossolo).  
- `geo_loc_name`, `lat_lon`, `depth`, `collection_date`, `elev`.  
- `target_gene = ITS`, `target_subfragment = ITS2`.  
- `pcr_primers = FWD:GCATCGATGAAGAACGCAGC;REV:TCCTCCGCTTATTGATATGC` (ITS3/ITS4).  
- `seq_meth = Illumina MiSeq`, `lib_layout = paired`.[file:6][file:5][web:24]

**A completar a partir de registros de campo e laboratório:**

- `lat_lon`, `collection_date`, `elev` específicos de cada amostra.  
- `nucl_acid_ext` (kit de extração utilizado, ex.: Qiagen DNeasy PowerSoil).  
- `pcr_cond` (ciclos e temperaturas de PCR usados no protocolo real).[file:4][file:7]

---

### 2.4. Ambiente de software (QIIME2, FastQC, MultiQC, R)

O pipeline foi executado em ambiente Linux (Ubuntu 22.04) utilizando Conda para gerenciar todas as dependências de bioinformática.[web:136][web:145]  
A instalação do QIIME2 seguiu o procedimento recomendado pela equipe do QIIME2, usando o ambiente “amplicon” oficial, que já inclui os plugins necessários (`q2-cutadapt`, `q2-dada2`, `q2-feature-classifier` etc.).[web:136][web:148]

Para a versão 2024.10, o ambiente foi criado com:

```bash
# baixar o YAML oficial da distribuição amplicon (Linux, py310)
wget https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-linux-conda.yml

# criar o ambiente conda
conda env create -n qiime2-amplicon-2024.10 \
  -f qiime2-amplicon-2024.10-py310-linux-conda.yml

# ativar
conda activate qiime2-amplicon-2024.10
```

Esse ambiente inclui automaticamente o núcleo do QIIME2 (`qiime2`) e os plugins padrão para análise de amplicon, incluindo `q2-cutadapt`, `q2-dada2` e `q2-feature-classifier`, evitando problemas de resolução de pacotes como `q2-cutadapt` em canais gerais do Conda.[web:136][web:145][web:148]

Ferramentas adicionais usadas no pipeline (FastQC, MultiQC, Cutadapt, VSEARCH) e o stack de R para análises downstream foram instalados sobre o mesmo ambiente QIIME2:

```bash
# ainda dentro do ambiente qiime2-amplicon-2024.10
conda install -c bioconda -c conda-forge \
  fastqc multiqc cutadapt vsearch

conda install -c conda-forge \
  r-base=4.3 r-tidyverse r-vegan r-phyloseq r-readr r-data.table r-ggplot2 r-patchwork r-optparse
```

Dessa forma, todas as etapas (QC com FastQC/MultiQC, processamento no QIIME2 e análises ecológicas em R) são executadas dentro de um único ambiente Conda reprodutível, baseado no YAML oficial da distribuição QIIME2 para amplicon.[web:136][web:140][web:145]
  

## 3. Controle de qualidade inicial (FastQC e MultiQC)

**Script:** `scripts/00-fastqc_multiqc.sh`  
**Pastas de saída:** `results/fastqc/`, `results/multiqc/`  

Antes da importação no QIIME2, todos os arquivos FASTQ brutos são avaliados com FastQC (vX.X) e os relatórios são sumarizados com MultiQC (vX.X).[web:66][web:67]  
São inspecionados: qualidade por base, distribuição de comprimento, conteúdo de GC, presença de adaptadores e sequências sobre‑representadas para cada amostra.

Script utilizado:

```bash
#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE}")/.." && pwd)"

mkdir -p "$PROJECT_ROOT/results/fastqc" "$PROJECT_ROOT/results/multiqc"

fastqc \
  "$PROJECT_ROOT"/data/raw/*/*.fastq* \
  --outdir "$PROJECT_ROOT/results/fastqc" \
  --threads 4

multiqc "$PROJECT_ROOT/results/fastqc" \
  --outdir "$PROJECT_ROOT/results/multiqc"
```

O relatório combinado `results/multiqc/multiqc_report.html` é utilizado para:

- Verificar se a qualidade por base se mantém alta ao longo dos 300 ciclos (espera‑se Phred ≥ 30 na maior parte do read).  
- Verificar se há contaminação significativa por adaptadores ou conteúdo excessivo de Ns.  
- Guiar a escolha de parâmetros de truncagem para DADA2.[file:7][cite:133]

O MultiQC indicou alertas esperados para bibliotecas de amplicon ITS2, principalmente em duplicação de sequências, sequências sobre-representadas e, em algumas amostras, conteúdo de adaptadores/primers. Esses sinais foram tratados nas etapas seguintes com Cutadapt e DADA2; a interpretação detalhada está em `README_results.qmd`.

---


## 4. Processamento bioinformático

### 4.1 Importação das leituras no QIIME 2

A análise das sequências ITS foi realizada no QIIME 2 (versão 2024.10, distribuição amplicon), utilizando o modo de importação via arquivo de manifest para leituras pareadas em formato FASTQ com codificação Phred33.[web:154][web:152] Os arquivos de leitura bruta foram organizados no diretório de trabalho do projeto em subpastas por sítio dentro de `data/raw/`, de forma que cada amostra possui dois arquivos (`*_1.fastq` e `*_2.fastq`) representando as leituras forward e reverse, respectivamente.[file:176]

```bash
cd ~/Documentos/1_Projetos/2026_UEFS_Projeto_metabarcoding_fungos_2026/02_Work/fungos_caatinga_V2

# organização dos FASTQ brutos
mkdir -p data/raw

# exemplo de estrutura esperada
ls data/raw
# CND  GBR  IBD  JGR  manifest.tsv

ls data/raw/CND
# CND_01_1.fastq  CND_01_2.fastq  CND_03_1.fastq  CND_03_2.fastq  ...
```

#### 4.1.1 Ajuste de caminhos do manifest e verificação

O conjunto original de amostras já possuía um arquivo de manifest compatível com o QIIME 2, contendo três colunas: `sample-id`, `forward-absolute-filepath` e `reverse-absolute-filepath`, no formato PairedEndFastqManifestPhred33V2.[file:176][web:166] Para reutilizar esse manifest no novo diretório de trabalho, o arquivo foi copiado de `data/processed/manifest.tsv` para `data/raw/manifest.tsv` e, em seguida, os caminhos absolutos internos foram atualizados com `sed`, substituindo o prefixo do projeto anterior pelo caminho do projeto atual.

```bash
cd ~/Documentos/1_Projetos/2026_UEFS_Projeto_metabarcoding_fungos_2026/02_Work/fungos_caatinga_V2

# copiar o manifest original para o local esperado pelo pipeline
cp data/processed/manifest.tsv data/raw/manifest.tsv

# atualizar os caminhos absolutos dos FASTQ no manifest
sed -i 's#/home/rmfsantos/Documentos/bioinfo/fungos-caatinga#/home/user/Documentos/1_Projetos/2026_UEFS_Projeto_metabarcoding_fungos_2026/02_Work/fungos_caatinga_V2#g' data/raw/manifest.tsv
```

Após a edição, a consistência do arquivo foi conferida visualmente, verificando-se se cada linha apontava para os arquivos FASTQ corretos dentro de `data/raw`.[file:176]

```bash
head -5 data/raw/manifest.tsv

# exemplo esperado de saída:
# sample-id    forward-absolute-filepath                                              reverse-absolute-filepath
# CND_01       /home/user/.../fungos_caatinga_V2/data/raw/CND/CND_01_1.fastq         /home/user/.../fungos_caatinga_V2/data/raw/CND/CND_01_2.fastq
```

Também foi verificada explicitamente a existência dos arquivos indicados na primeira linha do manifest para evitar erros de importação do tipo “Filepath on line X could not be found” no QIIME 2.[web:158][web:161]

```bash
ls data/raw/CND
ls data/raw/CND/CND_01_1.fastq
```

#### 4.1.2 Comando de importação no QIIME 2

Com o ambiente QIIME 2 ativado e o manifest validado, as leituras pareadas foram importadas para o formato interno `.qza` usando o tipo `SampleData[PairedEndSequencesWithQuality]` e o formato `PairedEndFastqManifestPhred33V2`, conforme recomendado na documentação oficial.[web:154][web:163]

```bash
conda activate qiime2-amplicon-2024.10

cd ~/Documentos/1_Projetos/2026_UEFS_Projeto_metabarcoding_fungos_2026/02_Work/fungos_caatinga_V2

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path data/raw/manifest.tsv \
  --output-path workflow/01-import/demux.qza \
  --input-format PairedEndFastqManifestPhred33V2
```

Em seguida, foi gerado um sumário de controle de qualidade das leituras importadas com `qiime demux summarize`, produzindo um arquivo `.qzv` explorável no QIIME 2 View.[web:152][web:157]

```bash
qiime demux summarize \
  --i-data workflow/01-import/demux.qza \
  --o-visualization results/qiime2/demux.qzv
```

Esse sumário foi utilizado posteriormente para inspecionar a distribuição do número de leituras por amostra e os perfis de qualidade Phred ao longo dos ciclos de sequenciamento, orientando a escolha dos parâmetros de truncamento e recorte na etapa de denoising com DADA2.[web:157]

O sumário `demux.qzv` indicou 50 amostras com leituras pareadas, totalizando 9.786.539 reads forward e 9.786.539 reads reverse antes de qualquer filtragem de qualidade.[web:194] O número de reads por amostra variou de 87.112 a 219.097, com mediana de 204.585 reads por amostra para ambas as direções, indicando profundidade de sequenciamento elevada e relativamente uniforme entre as amostras.[web:194] As distribuições de comprimento foram estreitas (medianas de 225 pb para leituras forward e 226 pb para reverse, interquartis ~223–228 pb), com medianas de qualidade Phred ≈ Q37–38 até ~200–210 pb e aumento gradual da dispersão apenas nos ciclos finais, padrão compatível com dados MiSeq de amplicons ITS.[web:194][web:195][web:204]

### 4.2 Remoção de primers ITS com Cutadapt (`scripts/01b-cutadapt.sh`)

Após a importação das leituras pareadas no QIIME 2, procedeu-se à remoção dos primers da região ITS utilizando o plugin `q2-cutadapt`, que encapsula o Cutadapt (https://cutadapt.readthedocs.io).[file:283] Essa etapa segue recomendações de boas práticas em metabarcoding fúngico para remover sequências de primers antes da inferência de ASVs, reduzindo vieses de abundância e melhorando a acurácia da anotação taxonômica.[file:283]

A remoção dos primers foi realizada pelo script `scripts/01b-cutadapt.sh`, executado a partir do diretório raiz do projeto. Esse script localiza automaticamente o diretório do projeto e aplica o corte de primers diretamente sobre o artefato desmultiplexado (`demux.qza`), gerando um novo conjunto de leituras pareadas aparadas (`demux-trimmed.qza`) e um sumário de qualidade correspondente (`demux-trimmed.qzv`):

```bash
#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE}")/.." && pwd)"

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences "$PROJECT_ROOT/workflow/01-import/demux.qza" \
  --p-front-f GCATCGATGAAGAACGCAGC \
  --p-front-r TCCTCCGCTTATTGATATGC \
  --p-error-rate 0.1 \
  --p-discard-untrimmed \
  --p-cores 0 \
  --o-trimmed-sequences "$PROJECT_ROOT/workflow/01b-cutadapt/demux-trimmed.qza" \
  --verbose

qiime demux summarize \
  --i-data "$PROJECT_ROOT/workflow/01b-cutadapt/demux-trimmed.qza" \
  --o-visualization "$PROJECT_ROOT/results/qiime2/demux-trimmed.qzv"
```

Forward e reverse primers ITS2 (`GCATCGATGAAGAACGCAGC` e `TCCTCCGCTTATTGATATGC`, respectivamente) foram especificados com `--p-front-f` e `--p-front-r`, admitindo taxa máxima de erro de 10% (`--p-error-rate 0.1`) e descartando leituras sem primers detectados em nenhuma das extremidades (`--p-discard-untrimmed`).[file:283] O parâmetro `--p-cores 0` permite o uso de todos os núcleos disponíveis para paralelização, tornando o processamento mais eficiente em máquinas com múltiplos cores.[file:283] As leituras aparadas resultantes (`demux-trimmed.qza`) foram utilizadas como entrada para inferência de ASVs com DADA2 e posterior classificação taxonômica com UNITE.[file:280][file:283]

### 4.3 Filtragem de qualidade e inferência de ASVs com DADA2

As leituras pareadas aparadas para remoção de primers ITS (`workflow/01b-cutadapt/demux-trimmed.qza`, Seção 4.2) foram submetidas à inferência de variantes de sequência de amplicon (ASVs) utilizando o plugin `q2-dada2` do QIIME 2.[file:283][web:332] O método `dada2 denoise-paired` executa filtragem de qualidade, modelagem de erros, desduplicação, junção de leituras pareadas e remoção de quimeras, produzindo uma tabela de frequências e um conjunto de sequências representativas de alta resolução.[file:283][web:332]

O processamento gerou a tabela de ASVs, as sequências representativas e as estatísticas de denoising no diretório `workflow/02-dada2`. O script correspondente ainda deve ser versionado em `scripts/`; o comando executado foi:

```bash
#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE}")/.." && pwd)"

mkdir -p "$PROJECT_ROOT/workflow/02-dada2"

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "$PROJECT_ROOT/workflow/01b-cutadapt/demux-trimmed.qza" \
  --p-trunc-len-f 220 \
  --p-trunc-len-r 220 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-max-ee-f 4 \
  --p-max-ee-r 4 \
  --p-trunc-q 0 \
  --p-n-threads 0 \
  --o-table "$PROJECT_ROOT/workflow/02-dada2/table.qza" \
  --o-representative-sequences "$PROJECT_ROOT/workflow/02-dada2/rep-seqs.qza" \
  --o-denoising-stats "$PROJECT_ROOT/workflow/02-dada2/denoising-stats.qza"
```

Os valores de truncamento (`trunc-len-f = 220` e `trunc-len-r = 220`) foram escolhidos com base na distribuição de comprimentos observada em `demux-trimmed.qzv`, na qual a maioria das leituras apresentou 219–228 nucleotídeos, garantindo remoção apenas da cauda de menor qualidade e uma ampla região de sobreposição entre as leituras forward e reverse.[web:320][web:331] O filtro de qualidade foi definido com número máximo de erros esperados igual a 4 em cada direção (`max-ee-f` e `max-ee-r`), um compromisso entre retenção de leituras e controle de erros, conforme recomendações para dados ITS com boa cobertura.[web:325][web:332] As ASVs obtidas (`table.qza` e `rep-seqs.qza`) foram utilizadas nas etapas subsequentes de atribuição taxonômica com o banco UNITE e análises de diversidade alfa e beta.[file:280][file:283]

#### 4.3.1 Geração de visualizações dos resultados do DADA2

Após a execução do DADA2, foram gerados sumários e visualizações dos resultados, incluindo: (i) tabela de abundância de ASVs por amostra, (ii) sequências representativas, e (iii) estatísticas detalhadas de denoising.[web:195][web:204]

```bash
# sumário da tabela de ASVs por amostra
qiime feature-table summarize \
  --i-table workflow/02-dada2/table.qza \
  --o-visualization results/qiime2/table.qzv \
  --m-sample-metadata-file data/external/metadata.tsv

# inspeção das sequências representativas
qiime feature-table tabulate-seqs \
  --i-data workflow/02-dada2/rep-seqs.qza \
  --o-visualization results/qiime2/rep-seqs.qzv

# estatísticas detalhadas do denoising
qiime metadata tabulate \
  --m-input-file workflow/02-dada2/denoising-stats.qza \
  --o-visualization results/qiime2/denoising-stats.qzv
```

Os arquivos `table.qzv` e `denoising-stats.qzv` foram utilizados para quantificar o número de ASVs, a distribuição de abundância de leituras por amostra após denoising e a fração de reads descartadas em cada etapa (filtragem de qualidade, inferência de erros, junção de pares e remoção de quimeras), informações que serão reportadas na seção de Resultados e no material suplementar do artigo.[web:195][file:169][web:200]

O sumário `denoising-stats.qzv` indicou que, após a filtragem de qualidade e correção de erros pelo DADA2, entre aproximadamente 72% e 94% das leituras por amostra passaram pelos filtros iniciais, com taxas de junção de pares geralmente entre 80% e 91% das leituras de entrada.[web:195][web:204] A fração final de leituras não quiméricas (non-chimeric) variou de cerca de 50% a 88% do total de reads de entrada, com a maioria das amostras apresentando retenção em torno de 70–80%, o que é consistente com valores reportados em estudos de metabarcoding fúngico em solos.[file:170][file:172][web:200]

### 4.4 Atribuição taxonômica com UNITE

As variantes de sequência de amplicon (ASVs) inferidas por DADA2 (`rep-seqs.qza`) foram classificadas taxonomicamente utilizando um classificador Naive Bayes pré-treinado no banco UNITE ITS para fungos, implementado no plugin `q2-feature-classifier` do QIIME 2.[file:283][web:351][web:354] Foi utilizado o classificador `unite_ver10_99_all_19.02.2025-Q2-2024.10.qza`, disponibilizado para a versão 2024.10 do QIIME 2 e armazenado em `data/external` no diretório do projeto.[web:355][web:364]

A atribuição taxonômica gravou as anotações em `workflow/03-taxonomy/taxonomy.qza`. O script correspondente ainda deve ser versionado em `scripts/`; o comando executado foi:

```bash
#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE}")/.." && pwd)"

mkdir -p "$PROJECT_ROOT/workflow/03-taxonomy" "$PROJECT_ROOT/results/qiime2"

qiime feature-classifier classify-sklearn \
  --i-classifier "$PROJECT_ROOT/data/external/unite_ver10_99_all_19.02.2025-Q2-2024.10.qza" \
  --i-reads "$PROJECT_ROOT/workflow/02-dada2/rep-seqs.qza" \
  --p-n-jobs 0 \
  --o-classification "$PROJECT_ROOT/workflow/03-taxonomy/taxonomy.qza"

qiime metadata tabulate \
  --m-input-file "$PROJECT_ROOT/workflow/03-taxonomy/taxonomy.qza" \
  --o-visualization "$PROJECT_ROOT/results/qiime2/taxonomy.qzv"
```

O método `classify-sklearn` aplica um classificador Naive Bayes pré-ajustado com limiar de confiança padrão 0,7 para limitar a profundidade taxonômica, conforme a documentação do QIIME 2 para dados ITS.[web:352][web:364] As anotações taxonômicas resultantes (`taxonomy.qza`) foram utilizadas, em conjunto com a tabela de ASVs, para gerar gráficos de barras de composição taxonômica e para agregação das abundâncias por diferentes níveis taxonômicos.[web:359][web:365]

### 4.5 Filtragem para ASVs fúngicas classificadas (`scripts/04-filtering.sh`)

Após a atribuição taxonômica, a tabela de ASVs e as sequências representativas foram filtradas para manter apenas variantes atribuídas ao reino Fungi e com classificação não indefinida (isto é, removendo rótulos "Unassigned").[web:360][web:362] Antes da filtragem taxonômica, a tabela e as sequências foram restringidas ao conjunto de ASVs presentes em `taxonomy.qza`, garantindo correspondência exata entre os identificadores de features e as entradas de taxonomia.[web:360][web:377]

Essas etapas foram automatizadas no script `scripts/04-filtering.sh`, que produz uma tabela final de ASVs fúngicas classificadas (`table_fungi_classified.qza`) e as sequências representativas correspondentes (`rep-seqs_fungi_classified.qza`):

```bash
#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE}")/.." && pwd)"

mkdir -p "$PROJECT_ROOT/workflow/04-filtering" "$PROJECT_ROOT/results/qiime2"

# 1) Garantir que apenas ASVs com taxonomia conhecida sejam mantidas
qiime feature-table filter-features \
  --i-table "$PROJECT_ROOT/workflow/02-dada2/table.qza" \
  --m-metadata-file "$PROJECT_ROOT/workflow/03-taxonomy/taxonomy.qza" \
  --o-filtered-table "$PROJECT_ROOT/workflow/04-filtering/table_taxa-matched.qza"

qiime feature-table filter-seqs \
  --i-data "$PROJECT_ROOT/workflow/02-dada2/rep-seqs.qza" \
  --m-metadata-file "$PROJECT_ROOT/workflow/03-taxonomy/taxonomy.qza" \
  --o-filtered-data "$PROJECT_ROOT/workflow/04-filtering/rep-seqs_taxa-matched.qza"

# 2) Filtrar a tabela e as sequências apenas para Fungi, removendo "Unassigned"
qiime taxa filter-table \
  --i-table "$PROJECT_ROOT/workflow/04-filtering/table_taxa-matched.qza" \
  --i-taxonomy "$PROJECT_ROOT/workflow/03-taxonomy/taxonomy.qza" \
  --p-include "Fungi" \
  --p-exclude "Unassigned" \
  --o-filtered-table "$PROJECT_ROOT/workflow/04-filtering/table_fungi_classified.qza"

qiime taxa filter-seqs \
  --i-sequences "$PROJECT_ROOT/workflow/04-filtering/rep-seqs_taxa-matched.qza" \
  --i-taxonomy "$PROJECT_ROOT/workflow/03-taxonomy/taxonomy.qza" \
  --p-include "Fungi" \
  --p-exclude "Unassigned" \
  --o-filtered-sequences "$PROJECT_ROOT/workflow/04-filtering/rep-seqs_fungi_classified.qza"

# 3) Sumário da tabela final de ASVs fúngicas
qiime feature-table summarize \
  --i-table "$PROJECT_ROOT/workflow/04-filtering/table_fungi_classified.qza" \
  --o-visualization "$PROJECT_ROOT/results/qiime2/table_fungi_classified.qzv"
```

A abordagem em dois passos (primeiro casar a tabela e as sequências com a taxonomia via `feature-table filter-features/seqs`, depois aplicar `taxa filter-table/seqs` para `include = "Fungi"` e `exclude = "Unassigned"`) segue as recomendações da documentação do QIIME 2 para evitar inconsistências de identificadores de features durante filtragens taxonômicas.[web:360][web:362] A tabela filtrada (`table_fungi_classified.qza`) foi utilizada como base para as etapas subsequentes de visualização da composição taxonômica (gráficos de barras) e para as análises de diversidade.[web:359][web:365]

# 4.2 Exportar a feature table filtrada e converter para TSV

# exporta a tabela em formato BIOM
mkdir -p data/processed/fungi/export
qiime tools export \
  --input-path workflow/04-filtering/table_fungi_classified.qza \
  --output-path data/processed/fungi/export

# converte o BIOM para TSV legível
biom convert \
  -i data/processed/fungi/export/feature-table.biom \
  -o data/processed/fungi/export/asv_table_fungi.tsv \
  --to-tsv

# O arquivo asv_table_fungi.tsv contém as ASVs (linhas) e as amostras (colunas),
# já filtradas apenas para sequências classificadas como k__Fungi.



### 4.6 Visualização da composição taxonômica (`scripts/05-taxa-barplot.sh`)

A composição taxonômica das comunidades fúngicas foi visualizada a partir da tabela filtrada de ASVs fúngicas classificadas (`table_fungi_classified.qza`) e das anotações taxonômicas obtidas com o classificador UNITE (`taxonomy.qza`).[web:359][web:365] Utilizou-se o plugin `q2-taxa` do QIIME 2 para gerar gráficos de barras interativos com a abundância relativa de táxons em diferentes níveis taxonômicos, integrando a tabela de metadados `data/external/metadata.tsv` para permitir a comparação entre supersítios e subáreas de amostragem.[file:286][web:362]

A geração dos barplots foi automatizada pelo script `scripts/05-taxa-barplot.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE}")/.." && pwd)"

mkdir -p "$PROJECT_ROOT/results/qiime2"

qiime taxa barplot \
  --i-table "$PROJECT_ROOT/workflow/04-filtering/table_fungi_classified.qza" \
  --i-taxonomy "$PROJECT_ROOT/workflow/03-taxonomy/taxonomy.qza" \
  --m-metadata-file "$PROJECT_ROOT/data/external/metadata.tsv" \
  --o-visualization "$PROJECT_ROOT/results/qiime2/taxa-barplot_fungi.qzv"
```

Os gráficos resultantes (`taxa-barplot_fungi.qzv`) permitiram inspecionar a distribuição relativa de filos, classes, ordens e gêneros fúngicos entre os diferentes supersítios da Caatinga, bem como explorar padrões de variação espacial na composição das comunidades.[web:359][web:365]

## Tabela de ASVs por supersítio com taxonomia

Após a filtragem para fungos e a classificação taxonômica, geramos uma tabela consolidada indicando, para cada ASV, em quais supersítios foi detectada (presença/ausência) e a respectiva anotação taxonômica.

1. Exportar a taxonomia e gerar `taxonomy_fungi.tsv` (todas as ASVs ITS anotadas como fungos):

```bash
mkdir -p data/processed/fungi/taxonomy_export

qiime tools export \
  --input-path workflow/03-taxonomy/taxonomy.qza \
  --output-path data/processed/fungi/taxonomy_export

python3 - << 'PY'
import pandas as pd
from pathlib import Path

inp = Path("data/processed/fungi/taxonomy_export/taxonomy.tsv")
out = Path("data/processed/fungi/taxonomy_fungi.tsv")

df = pd.read_csv(inp, sep="\t")
# renomeia 'Feature ID' -> 'Feature.ID' para compatibilizar com outros scripts
if "Feature ID" in df.columns:
    df = df.rename(columns={"Feature ID": "Feature.ID"})
df.to_csv(out, sep="\t", index=False)
print(f"taxonomy_fungi.tsv salvo em: {out}, n_linhas =", len(df))
PY
```

2. A partir das listas `ASVs_<Supersitio>.txt` (geradas previamente a partir da tabela `table_fungi_classified.qza`), construímos a matriz presença/ausência ASV × supersítio:

```bash
bash scripts/monta_tabela_asv_supersitio.sh
# saída: results/asv_por_supersitio.tsv
# colunas: ASV, Canudos, Gruta_Brejoes, Libraba, Jaguara, ...
```

3. Em seguida, anexamos a taxonomia de cada ASV usando `taxonomy_fungi.tsv`:

```bash
bash scripts/anexa_taxonomia_asv_supersitio.sh
# saída final: results/asv_por_supersitio_com_taxonomia.tsv
# colunas: ASV, supersítios (1 = presente, 0 = ausente) e Taxon
```

O arquivo `results/asv_por_supersitio_com_taxonomia.tsv` é utilizado pela doutoranda para explorar ASVs exclusivas ou compartilhadas entre supersítios, bem como para sumarizar padrões taxonômicos (por filo, gênero, etc.).

## Etapa 06 — Filogenia exploratória

A etapa 06 do pipeline realizou a inferência filogenética exploratória das ASVs fúngicas no QIIME 2 a partir do arquivo `workflow/04-filtering/rep-seqs_fungi.qza`. Foi utilizado o comando `qiime phylogeny align-to-tree-mafft-fasttree`, que executa o alinhamento múltiplo com MAFFT, o mascaramento do alinhamento, a inferência da árvore com FastTree e o enraizamento pelo método do ponto médio (midpoint rooting) [cite:8][cite:9][cite:15].

Os artefatos gerados foram salvos em `workflow/05-phylogeny/` e copiados para `results/qiime2/`, incluindo `aligned-rep-seqs_fungi.qza`, `masked-aligned-rep-seqs_fungi.qza`, `unrooted-tree_fungi.qza` e `rooted-tree_fungi.qza`. Além disso, as árvores foram exportadas com `qiime tools export` para o formato Newick (`tree.nwk`) em `results/qiime2/exported-tree/`, permitindo visualização externa em programas como iTOL, FigTree ou análise posterior no R [cite:27][cite:47][cite:49].

Como o marcador utilizado foi ITS2, a árvore foi tratada como uma análise **exploratória**, útil para inspeção geral do agrupamento das ASVs e apoio a análises complementares, mas com interpretação cautelosa. Isso porque o ITS de fungos pode ser informativo para grupos próximos, porém tende a não produzir alinhamentos e filogenias robustos quando o conjunto inclui táxons muito divergentes [cite:40][cite:6].

### Script executado

```bash
bash scripts/06-phylogeny.sh
```

### Comando principal da etapa

```bash
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences workflow/04-filtering/rep-seqs_fungi.qza \
  --o-alignment workflow/05-phylogeny/aligned-rep-seqs_fungi.qza \
  --o-masked-alignment workflow/05-phylogeny/masked-aligned-rep-seqs_fungi.qza \
  --o-tree workflow/05-phylogeny/unrooted-tree_fungi.qza \
  --o-rooted-tree workflow/05-phylogeny/rooted-tree_fungi.qza \
  --p-n-threads auto
```

### Exportação das árvores

```bash
qiime tools export \
  --input-path workflow/05-phylogeny/unrooted-tree_fungi.qza \
  --output-path results/qiime2/exported-tree/unrooted-tree

qiime tools export \
  --input-path workflow/05-phylogeny/rooted-tree_fungi.qza \
  --output-path results/qiime2/exported-tree/rooted-tree
```

## 5. Geração de figuras para o manuscrito

As figuras do manuscrito foram planejadas para serem geradas a partir dos artefatos finais do pipeline em QIIME 2 e das tabelas derivadas em R, com foco em reprodutibilidade, padronização visual e exportação em formatos adequados para submissão editorial.[web:27][web:189] Os arquivos finais foram organizados no diretório `results/figures/`, com versões em PNG e TIFF em alta resolução (600 dpi), apropriadas para uso em manuscritos, apresentações e material suplementar.[web:186][web:189]

A preparação das figuras seguiu duas estratégias complementares. Primeiro, os artefatos `.qza` e `.qzv` produzidos no QIIME 2 foram exportados com `qiime tools export` sempre que necessário para obtenção de arquivos tabulares intermediários utilizáveis em R.[web:27] Segundo, os gráficos finais foram montados com `ggplot2`, permitindo controle explícito de cores, rótulos, escalas e dimensões de exportação.[web:186][web:189]

### 5.1. Figura 3 — Qualidade e retenção de reads após DADA2

A Figura 3 foi definida como um painel voltado à documentação do desempenho da etapa de denoising, com dois componentes: (A) distribuição do número de leituras por amostra antes e após o DADA2; e (B) porcentagem de retenção de reads nas etapas de filtragem, denoising, junção de pares e remoção de quimeras.[web:155][web:153]

Os dados para essa figura foram extraídos a partir do artefato `workflow/02-dada2/denoising-stats.qza`, exportado com `qiime tools export`, o que gera um arquivo tabular `stats.tsv` contendo, para cada amostra, o número de reads de entrada (`input`), reads filtradas (`filtered`), reads denoised, reads mescladas (`merged`) e reads não-quiméricas (`non-chimeric`). Esse procedimento é compatível com o fluxo recomendado na documentação do QIIME 2 para reutilização dos dados em análises externas e visualizações personalizadas em R.[web:27][web:155][web:153]

Quando necessário, a tabela de frequências por amostra após o DADA2 também pode ser obtida a partir de `workflow/02-dada2/table.qza` usando `qiime feature-table summarize` ou `qiime feature-table summarize-plus`, sendo este último capaz de produzir frequências por amostra e por feature em formato de metadados exportáveis.[web:157][web:160]

#### Exportação dos dados para a Figura 3

```bash
mkdir -p data/processed/qc_export

qiime tools export \
  --input-path workflow/02-dada2/denoising-stats.qza \
  --output-path data/processed/qc_export/denoising_stats_export
```

Em versões do QIIME 2 que disponham do plugin `feature-table summarize-plus`, a extração de frequências por amostra pode ser feita adicionalmente com:

```bash
qiime feature-table summarize-plus \
  --i-table workflow/02-dada2/table.qza \
  --o-feature-frequencies data/processed/qc_export/feature-freq.qza \
  --o-sample-frequencies data/processed/qc_export/sample-freq.qza \
  --o-summary data/processed/qc_export/table_summary.qzv
```

#### Geração dos gráficos em R

A montagem final da Figura 3 foi realizada em R com `ggplot2`, a partir do arquivo `data/processed/qc_export/denoising_stats_export/stats.tsv` e dos metadados amostrais (`data/external/metadata.tsv`). O script correspondente foi salvo em `scripts/fig3_qc_dada2.R`.[web:182][web:186]

A Figura 3A corresponde a um boxplot do número de reads por amostra antes do DADA2 (`input`) e após a remoção de quimeras (`non-chimeric`). A Figura 3B corresponde a boxplots da porcentagem de retenção em cada etapa do processamento, calculada em relação ao total de reads de entrada por amostra.[web:155][web:153]

A exportação dos gráficos foi feita em PNG e TIFF com 600 dpi, resolução comumente aceita para publicação científica.[web:189]

```r
ggsave(
  filename = file.path(fig_dir, "Fig3_QC_DADA2_reads_retention.png"),
  plot = fig3,
  width = 16, height = 8, units = "cm", dpi = 600
)

ggsave(
  filename = file.path(fig_dir, "Fig3_QC_DADA2_reads_retention.tiff"),
  plot = fig3,
  width = 16, height = 8, units = "cm", dpi = 600,
  compression = "lzw"
)
```

#### Script executado

```bash
Rscript scripts/fig3_qc_dada2.R
```

#### Saídas esperadas

```text
results/figures/Fig3_QC_DADA2_reads_retention.png
results/figures/Fig3_QC_DADA2_reads_retention.tiff
```

### 5.2. Padronização para as demais figuras

As demais figuras do manuscrito seguirão a mesma lógica de organização: (i) exportação dos dados relevantes dos artefatos QIIME 2; (ii) preparação de tabelas derivadas em `data/processed/`; (iii) geração dos gráficos em R com scripts versionados em `scripts/`; e (iv) salvamento final em `results/figures/` com dimensões e resolução adequadas para publicação.[web:27][web:157][web:189]

Esse padrão será aplicado às próximas figuras do artigo, incluindo composição taxonômica por supersítio, diversidade alfa e beta, filogenia exploratória e compartilhamento de ASVs entre áreas de amostragem.[web:182][web:189]

## Geração da Figura 4: composição taxonômica

A Figura 4 foi gerada com um script específico (`fig4`), desenvolvido para resumir e visualizar a composição taxonômica das comunidades fúngicas a partir da tabela de abundância e da classificação taxonômica previamente obtidas no pipeline de bioinformática. O objetivo do script foi produzir gráficos de barras empilhadas (*stacked bar plots*) mostrando a abundância relativa dos principais grupos taxonômicos em diferentes níveis hierárquicos, com foco em **filo** (Fig. 4A) e **classe** (Fig. 4B), uma estratégia amplamente empregada para visualização da composição global de comunidades fúngicas em estudos de metabarcoding [file:6][file:241].

O script `fig4` utilizou como entrada uma tabela de abundância por amostra e a tabela de taxonomia associada às ASVs fúngicas. A partir desses arquivos, as abundâncias foram agregadas por nível taxonômico, convertidas em abundância relativa e organizadas em formato apropriado para plotagem. Em seguida, o script gerou gráficos de barras empilhadas com cores distintas para cada táxon, permitindo comparar visualmente a contribuição relativa dos grupos dominantes entre as amostras ou grupos experimentais [file:239][file:241].

Na etapa de composição visual, o script foi configurado para destacar apenas os táxons mais representativos em cada nível taxonômico, agrupando os demais em uma categoria residual, como `Outros`, para reduzir ruído visual e facilitar a interpretação. Essa abordagem é recomendada quando há grande número de táxons raros, condição comum em dados de metabarcoding, e torna a figura mais adequada para apresentação em README, relatórios e manuscritos [file:239][file:6].

Os produtos finais do script corresponderam às imagens da Figura 4A e Figura 4B, representando, respectivamente, a composição taxonômica em nível de filo e classe. Essas figuras foram exportadas em formato de imagem para inserção direta no README e servem como resumo visual da estrutura taxonômica das comunidades fúngicas recuperadas após o processamento e a classificação das ASVs [file:245][file:246].

### Finalidade do script `fig4`

- Agregar abundâncias de ASVs por nível taxonômico.
- Converter contagens em abundância relativa por amostra ou tratamento.
- Selecionar os grupos taxonômicos mais abundantes.
- Gerar gráficos de barras empilhadas para composição taxonômica.
- Exportar as figuras finais para uso no README e em materiais de divulgação dos resultados.

### Saídas esperadas

- `Fig4A_taxonomic_composition_phylum.jpg`
- `Fig4B_taxonomic_composition_class.jpg`

### Geração das figuras de composição taxonômica

As figuras de composição taxonômica foram geradas a partir da tabela de abundância de ASVs e da tabela de classificação taxonômica, após filtragem para Fungi e integração com os metadados amostrais. Inicialmente, os dados foram processados no QIIME 2 para inspeção da composição relativa por amostra com `qiime taxa barplot`, permitindo verificar a distribuição dos principais táxons recuperados ao longo das amostras e superstítios [file:239][file:240][file:6].

Para a elaboração das figuras finais utilizadas no manuscrito, os dados exportados do QIIME 2 foram importados no R e organizados em objetos `phyloseq`, o que permitiu agregar as abundâncias em diferentes níveis taxonômicos e produzir gráficos com `ggplot2` [file:239][file:240]. A composição taxonômica foi resumida principalmente nos níveis de filo e gênero, expressa em abundância relativa média por superstítio, abordagem amplamente utilizada em estudos de comunidades fúngicas do solo [file:239][file:240].

No nível de gênero, os táxons foram agregados, transformados em abundância relativa e resumidos por superstítio. Para melhorar a legibilidade da figura, apenas os gêneros mais abundantes foram mantidos individualmente, enquanto os demais foram agrupados na categoria `Others`; registros sem classificação em nível de gênero foram reunidos em `Unclassified` [file:239][file:240]. No nível de filo, foi aplicado o mesmo princípio, permitindo comparar padrões mais amplos de composição taxonômica entre as áreas estudadas [file:239][file:240].

A figura final de composição taxonômica foi construída como uma **figura combinada**, reunindo em um único painel o barplot de abundância relativa em nível de filo e o barplot em nível de gênero. Essa montagem foi realizada no R com o pacote `patchwork`, que permite combinar múltiplos gráficos produzidos em `ggplot2` em um arranjo único, padronizado e apropriado para publicação [file:239]. Os gráficos foram exportados em alta resolução nos formatos PNG e PDF.

#### Comandos no QIIME 2

```bash
qiime taxa barplot \
  --i-table workflow/04-filtering/table_fungi_classified.qza \
  --i-taxonomy results/qiime2/taxonomy.qza \
  --m-metadata-file data/processed/sample-metadata.tsv \
  --o-visualization results/qiime2/taxa-barplot_fungi.qzv
```

```bash
qiime metadata tabulate \
  --m-input-file results/qiime2/taxonomy.qza \
  --o-visualization results/qiime2/taxonomy.qzv
```

```bash
qiime tools export \
  --input-path workflow/04-filtering/table_fungi_classified.qza \
  --output-path results/tables/exported-table_fungi_classified
```

```bash
biom convert \
  -i results/tables/exported-table_fungi_classified/feature-table.biom \
  -o results/tables/exported-table_fungi_classified/feature-table.tsv \
  --to-tsv
```

#### Script da figura combinada no R

```r
library(phyloseq)
library(ggplot2)
library(dplyr)
library(forcats)
library(scales)
library(RColorBrewer)
library(patchwork)

# objeto phyloseq já filtrado para fungos
# ps_fungi

group_var <- "Supersite"
site_order <- c("Canudos", "Gruta_Brejoes", "Ibiraba", "Jaguara")

# ----------------------------
# FILO
# ----------------------------
ps_phylum <- tax_glom(ps_fungi, taxrank = "Phylum", NArm = FALSE)
tax_table(ps_phylum)[, "Phylum"] <- as.character(tax_table(ps_phylum)[, "Phylum"])
tax_table(ps_phylum)[is.na(tax_table(ps_phylum)[, "Phylum"]), "Phylum"] <- "Unclassified"
tax_table(ps_phylum)[tax_table(ps_phylum)[, "Phylum"] == "", "Phylum"] <- "Unclassified"

ps_phylum_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x))
df_phylum <- psmelt(ps_phylum_rel)

df_phylum_site <- df_phylum %>%
  group_by(.data[[group_var]], Phylum) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

top_phyla <- df_phylum_site %>%
  group_by(Phylum) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice_head(n = 10) %>%
  pull(Phylum)

df_phylum_site <- df_phylum_site %>%
  mutate(Phylum_plot = ifelse(Phylum %in% top_phyla, Phylum, "Others")) %>%
  group_by(.data[[group_var]], Phylum_plot) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

df_phylum_site[[group_var]] <- factor(df_phylum_site[[group_var]], levels = site_order)

phylum_order <- df_phylum_site %>%
  group_by(Phylum_plot) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Phylum_plot)

df_phylum_site$Phylum_plot <- factor(df_phylum_site$Phylum_plot, levels = phylum_order)

phylum_cols <- c(
  "Ascomycota" = "#E69F00",
  "Basidiomycota" = "#56B4E9",
  "Mortierellomycota" = "#009E73",
  "Mucoromycota" = "#F0E442",
  "Glomeromycota" = "#0072B2",
  "Chytridiomycota" = "#D55E00",
  "Rozellomycota" = "#CC79A7",
  "Unclassified" = "plum3",
  "Others" = "grey70"
)

missing_phylum_cols <- setdiff(levels(df_phylum_site$Phylum_plot), names(phylum_cols))
if (length(missing_phylum_cols) > 0) {
  extra_cols <- colorRampPalette(brewer.pal(8, "Set3"))(length(missing_phylum_cols))
  names(extra_cols) <- missing_phylum_cols
  phylum_cols <- c(phylum_cols, extra_cols)
}

p_phylum <- ggplot(df_phylum_site,
                   aes(x = .data[[group_var]], y = Abundance, fill = Phylum_plot)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = phylum_cols) +
  labs(x = NULL, y = "Relative abundance (%)", fill = "Phylum") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "right"
  )

# ----------------------------
# GÊNERO
# ----------------------------
ps_genus <- tax_glom(ps_fungi, taxrank = "Genus", NArm = FALSE)
tax_table(ps_genus)[, "Genus"] <- as.character(tax_table(ps_genus)[, "Genus"])
tax_table(ps_genus)[is.na(tax_table(ps_genus)[, "Genus"]), "Genus"] <- "Unclassified"
tax_table(ps_genus)[tax_table(ps_genus)[, "Genus"] == "", "Genus"] <- "Unclassified"

ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))
df_genus <- psmelt(ps_genus_rel)

df_genus_site <- df_genus %>%
  group_by(.data[[group_var]], Genus) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

top_genera <- df_genus_site %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  slice_head(n = 15) %>%
  pull(Genus)

df_genus_site <- df_genus_site %>%
  mutate(Genus_plot = ifelse(Genus %in% top_genera, Genus, "Others")) %>%
  group_by(.data[[group_var]], Genus_plot) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

df_genus_site[[group_var]] <- factor(df_genus_site[[group_var]], levels = site_order)

genus_order <- df_genus_site %>%
  group_by(Genus_plot) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Genus_plot)

genus_order <- c(setdiff(genus_order, c("Others", "Unclassified")),
                 intersect(c("Others", "Unclassified"), genus_order))

df_genus_site$Genus_plot <- factor(df_genus_site$Genus_plot, levels = genus_order)

genus_cols <- colorRampPalette(c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
  "#66a61e", "#e6ab02", "#a6761d", "#666666"
))(length(levels(df_genus_site$Genus_plot)))
names(genus_cols) <- levels(df_genus_site$Genus_plot)

if ("Others" %in% names(genus_cols)) genus_cols["Others"] <- "grey70"
if ("Unclassified" %in% names(genus_cols)) genus_cols["Unclassified"] <- "plum3"

p_genus <- ggplot(df_genus_site,
                  aes(x = .data[[group_var]], y = Abundance, fill = Genus_plot)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = genus_cols) +
  labs(x = NULL, y = "Relative abundance (%)", fill = "Genus") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    legend.position = "right"
  )

# ----------------------------
# FIGURA COMBINADA
# ----------------------------
p_combined <- p_phylum / p_genus +
  plot_annotation(tag_levels = "A")

ggsave("Fig_taxonomic_composition_combined.png",
       p_combined, width = 12, height = 14, dpi = 300, bg = "white")

ggsave("Fig_taxonomic_composition_combined.pdf",
       p_combined, width = 12, height = 14, bg = "white")

p_combined
```
Legenda sugerida: Figura X. Composição taxonômica das comunidades fúngicas do solo nos superstítios Canudos, Gruta_Brejoes, Ibiraba e Jaguara, expressa como abundância relativa média. O painel A mostra a composição em nível de filo e o painel B a composição em nível de gênero, com agrupamento dos táxons menos abundantes em Others e dos registros sem classificação em Unclassified
A abundância relativa média foi então calculada por superstítio, utilizando a variável de agrupamento presente nos metadados amostrais. Para facilitar a interpretação visual, foram mantidos no gráfico apenas os 15 gêneros mais abundantes no conjunto total de amostras, enquanto os demais foram reunidos na categoria `Others`; registros sem identificação em nível de gênero foram mantidos como `Unclassified` [file:240][file:239].

## Controle de qualidade do processamento no DADA2

O processamento das sequências pareadas foi realizado no QIIME 2 utilizando o plugin `q2-dada2` (`qiime dada2 denoise-paired`), que executa simultaneamente os passos de filtragem de qualidade, correção de erros, junção de leituras pareadas e remoção de quimeras. A saída principal desse passo inclui a tabela de abundância de ASVs (`table.qza`), as sequências representativas (`rep-seqs.qza`) e as estatísticas de desnoisificação por amostra (`denoising-stats.qza`). [web:385][web:454]

Para resumir o desempenho do DADA2, o artefato `denoising-stats.qza` foi exportado para um arquivo tabular (`dada2-stats.tsv`) usando o comando `qiime tools export`, seguindo a documentação oficial do QIIME 2. Esse arquivo contém, para cada amostra, o número de leituras em cada etapa (`input`, `filtered`, `denoised`, `merged` e `non-chimeric`) e as porcentagens correspondentes em relação ao número de leituras de entrada. Essas estatísticas são amplamente utilizadas para interpretar a eficiência da etapa de desnoisificação e da remoção de quimeras. [web:27][web:381][web:383]

A partir de `dada2-stats.tsv`, desenvolvemos um script em R (`scripts/figure_combined_dada2.R`) para gerar uma figura de controle de qualidade em dois painéis. No painel (A), os números de leituras de cada amostra foram organizados em formato longo e representados como linhas que conectam as etapas sucessivas do DADA2 (`input`, `filtered`, `denoised`, `merged` e `non-chimeric`). Esse painel permite visualizar a trajetória de perda de leituras ao longo do pipeline, identificando amostras que sofrem quedas acentuadas em etapas específicas (por exemplo, filtragem, junção das leituras ou remoção de quimeras). A escala vertical foi formatada com os rótulos abreviados (K para milhares) utilizando a função `label_number(scale_cut = cut_short_scale())` do pacote `scales`. [web:403][web:409][web:452]

No painel (B), utilizamos a coluna `percentage of input non-chimeric` para calcular a porcentagem de leituras de entrada que permaneceram como não quiméricas ao final do processamento, ordenando as amostras do menor para o maior valor e representando essa métrica em um gráfico de barras horizontais. Essa abordagem, frequentemente recomendada em tutoriais e discussões sobre interpretação do `denoising-stats` do DADA2, fornece uma visão global da retenção final de leituras por amostra e facilita a identificação de bibliotecas com desempenho atípico (baixa retenção), que podem demandar reavaliação de parâmetros de truncamento ou filtragem. [web:381][web:383][web:341]

Os gráficos foram produzidos com o pacote `ggplot2` e combinados em uma única figura etiquetada como painéis A e B por meio do pacote `patchwork`, seguindo uma prática comum para arranjos de múltiplos painéis em estudos de microbioma. A figura resultante foi exportada em alta resolução nos formatos `.png` e `.pdf` e utilizada como sumário visual do controle de qualidade do processamento no DADA2, antes das análises de composição taxonômica e diversidade. [web:452][web:345]

**Figura X.** Controle de qualidade do processamento no DADA2. (A) Contagem de leituras por amostra em cada etapa do DADA2 (`input`, `filtered`, `denoised`, `merged` e `non-chimeric`). (B) Porcentagem de leituras de entrada retidas como leituras não quiméricas ao final do processamento para cada amostra.

Figure 3. Quality-control summary of DADA2 processing.
(A) Read counts per sample across the main DADA2 processing steps (input, filtered, denoised, merged, and non-chimeric). This panel summarizes read loss throughout quality filtering, denoising, paired-read merging, and chimera removal.
(B) Percentage of input reads retained as non-chimeric reads after DADA2 for each sample, ordered from lowest to highest retention. This panel highlights variation in final read retention among samples and allows identification of libraries with comparatively low recovery after denoising.

Versão mais curta
Figure 3. DADA2 quality-control summary.
(A) Read counts per sample across DADA2 processing steps.
(B) Percentage of input reads retained as non-chimeric reads after DADA2, ordered by sample.

## Scripts de processamento e reprodutibilidade

Todo o processamento bioinformático foi organizado em scripts versionados no diretório `workflow/`, de forma a permitir a reexecução completa do pipeline a partir dos arquivos de leitura bruta. As etapas principais incluíram: (i) importação das leituras para o QIIME 2 (`01-import.sh`), (ii) remoção de adaptadores/trimming inicial quando aplicável, (iii) desnoisificação com DADA2 (`02-dada2`), (iv) atribuição taxonômica, e (v) geração de tabelas e figuras de composição taxonômica. [web:385][web:345]

No contexto da etapa de desnoisificação, o script `workflow/02-dada2` produz o artefato `denoising-stats.qza`, a partir do qual as estatísticas de qualidade por amostra são exportadas para o arquivo tabular `results/qiime2/dada2-stats.tsv` utilizando o comando `qiime tools export`, conforme recomendado na documentação do QIIME 2. Esse arquivo consolida, para cada amostra, as contagens de leituras e as porcentagens associadas em cada uma das etapas do DADA2. [web:27][web:381]

Os gráficos de controle de qualidade do DADA2 apresentados na Figura X foram gerados com o script em R `scripts/figure_combined_dada2.R`, que lê o arquivo `results/qiime2/dada2-stats.tsv`, realiza a limpeza e conversão dos campos numéricos, e constrói dois painéis com `ggplot2` e `patchwork`: (A) a trajetória de contagens de leituras por amostra nas etapas `input`, `filtered`, `denoised`, `merged` e `non-chimeric`, e (B) a distribuição da porcentagem de leituras de entrada retidas como não quiméricas para cada amostra. A figura final é salva automaticamente em `results/figures/` nos formatos `.png` e `.pdf`, garantindo a rastreabilidade entre os artefatos do QIIME 2, os scripts de análise em R e as figuras utilizadas no manuscrito. [web:403][web:409][web:452]

## Atualizações de hoje

Hoje foi estruturado o fluxo de integração entre a tabela de abundância de ASVs fúngicos, a taxonomia expandida e o metadata das amostras, com foco em gerar tabelas prontas para análises por amostra e por supersítio.[cite:614][cite:571]

### 1. Alinhamento entre abundância e taxonomia

Partiu-se de dois objetos principais no R: `asv_tab_fungi` (matriz de abundância de ASVs fúngicos) e `tax_tidy` (taxonomia expandida por ASV). A interseção entre `rownames(asv_tab_fungi)` e `tax_tidy$ASV_ID` foi usada para manter apenas ASVs presentes em ambos os objetos, garantindo compatibilidade entre abundância e anotação taxonômica.[cite:614][cite:617]

Em seguida, os objetos foram reordenados com base no vetor `common_asvs`, e a equivalência entre `tax_tidy2$ASV_ID` e `rownames(asv_tab_filt)` foi validada com `stopifnot(identical(...))`. Esse passo assegura que cada linha da matriz de abundância corresponde exatamente à mesma ASV na tabela taxonômica.[cite:614]

### 2. Geração da tabela consolidada ASV + taxonomia

Após o alinhamento, foi criada a tabela `asv_tax_table_fungi`, combinando as colunas taxonômicas com a matriz de abundância por amostra. Esse tipo de tabela preserva a resolução por ASV e ao mesmo tempo facilita sumarizações por níveis taxonômicos, como phylum, class ou genus.[cite:616][cite:622]

O arquivo exportado foi:

```r
results/asv_table_fungi_with_taxonomy.tsv
```

Esse arquivo deve ser usado como tabela-base para análises por amostra, incluindo riqueza, diversidade alfa, diversidade beta, ordenações e testes multivariados em função do metadata.[cite:617][cite:623]

### 3. Leitura e uso do metadata

Foi identificado que o arquivo de metadata está na raiz do projeto com o nome `metadata.tsv`, em conformidade com o formato tabular esperado pelo ecossistema QIIME 2 e por fluxos de análise em R.[cite:571][cite:564]

A leitura foi feita com `read.delim(..., row.names = 1)`, de modo que os IDs das amostras no metadata possam ser comparados diretamente com os nomes das colunas da matriz de abundância. Esse passo é essencial para agrupar amostras por supersítio e para análises dependentes de variáveis experimentais ou ambientais.[cite:571][cite:579]

### 4. Agregação da abundância por supersítio

Com base na coluna de supersítio presente no `metadata.tsv`, a matriz ASV x amostra foi convertida para formato longo, associada ao metadata, e agregada por `ASV_ID` e `Supersitio` usando soma das abundâncias. O resultado foi uma nova matriz ASV x supersítio, em que cada coluna representa a abundância total da ASV em um supersítio.[cite:619][cite:620]

O arquivo exportado foi:

```r
results/asv_por_supersitio_com_taxonomia.tsv
```

Esse arquivo é indicado para descrições agregadas da comunidade fúngica por supersítio, como composição taxonômica total, riqueza acumulada, distribuição de grupos dominantes e comparação geral entre supersítios.[cite:620][cite:617]

### 5. Conversão para presença/ausência

A tabela agregada por supersítio também foi convertida para formato binário, atribuindo 1 para presença e 0 para ausência de cada ASV em cada supersítio. Isso gerou uma segunda versão voltada para análises baseadas em ocorrência, em vez de abundância.[cite:618][cite:626]

O arquivo exportado foi:

```r
results/asv_por_supersitio_presenca_ausencia.tsv
```

Essa tabela é útil para quantificar ASVs exclusivos e compartilhados entre supersítios, estimar riqueza observada e calcular distâncias baseadas em presença/ausência, como Jaccard.[cite:618][cite:626]

### 6. Interpretação dos três arquivos gerados

Os arquivos gerados hoje têm funções complementares no pipeline analítico:[cite:617][cite:623]

- `asv_table_fungi_with_taxonomy.tsv`: tabela principal por amostra, adequada para análises ecológicas clássicas com abundância por amostra.[cite:617][cite:623]
- `asv_por_supersitio_com_taxonomia.tsv`: tabela agregada por supersítio, adequada para sínteses de composição e comparação entre supersítios.[cite:620]
- `asv_por_supersitio_presenca_ausencia.tsv`: tabela binária por supersítio, adequada para análises de compartilhamento, exclusividade e métricas baseadas em ocorrência.[cite:618][cite:626]

### 7. Preparação do gráfico de barras por phylum

Também foi montado o script para gerar um gráfico de barras empilhadas com abundância relativa por phylum em cada supersítio, usando a tabela `asv_por_supersitio_com_taxonomia.tsv`. Para isso, as abundâncias foram somadas por combinação de `Supersitio` e `Phylum`, e depois transformadas em abundância relativa dentro de cada supersítio.[cite:631][cite:632]

Esse tipo de visualização permite comparar a composição taxonômica entre supersítios em escala percentual, que é a forma mais comum de exibir dados composicionais em metabarcoding.[cite:632][cite:638]

### 8. Próximos usos recomendados

A partir dos arquivos produzidos hoje, o fluxo pode seguir em três direções principais:[cite:617][cite:623]

- análises por amostra: riqueza, Shannon, Bray-Curtis, NMDS/PCoA e PERMANOVA usando `asv_table_fungi_with_taxonomy.tsv` + `metadata.tsv`.[cite:623][cite:626]
- análises por supersítio: composição por phylum/class, riqueza total e ASVs exclusivas usando `asv_por_supersitio_com_taxonomia.tsv` e `asv_por_supersitio_presenca_ausencia.tsv`.[cite:620][cite:618]
- produção de figuras para resultados: barras empilhadas de abundância relativa, tabelas-resumo e gráficos comparativos entre supersítios.[cite:631][cite:632]

## Análise de diversidade alfa

A diversidade alfa foi calculada em nível de amostra a partir da tabela de abundância de ASVs fúngicos com taxonomia (`results/asv_table_fungi_with_taxonomy.tsv`) em conjunto com o arquivo de metadados das amostras (`metadata.tsv`). O arquivo de metadata foi lido como uma tabela tabulada por tabulação, com os identificadores das amostras armazenados como nomes de linha, em conformidade com o formato comumente utilizado em fluxos de trabalho do QIIME 2.[cite:563][cite:574]

### Tabelas de entrada e estrutura dos dados

A tabela `results/asv_table_fungi_with_taxonomy.tsv` contém colunas de anotação taxonômica seguidas pelas colunas de abundância de cada amostra individual. As porções taxonômica e de abundância foram separadas identificando-se a coluna `Confidence`, de modo que apenas a matriz de contagens foi utilizada nos cálculos de diversidade.[cite:564]

Em seguida, a matriz de abundância foi transposta para a orientação padrão de matriz comunidade exigida pelo pacote `vegan`, com amostras nas linhas e ASVs nas colunas.[cite:725][cite:726]

### Métricas de diversidade alfa

A riqueza observada foi calculada com `vegan::specnumber`, que retorna o número de táxons observados por amostra.[cite:725][cite:726] O índice de Shannon foi calculado com `vegan::diversity(..., index = "shannon")`, resumindo simultaneamente componentes de riqueza e equabilidade dentro de cada comunidade amostral.[cite:727][cite:729] O índice de Simpson também foi calculado com `vegan::diversity(..., index = "simpson")` como uma métrica adicional de diversidade intra-amostra.[cite:726]

Os índices resultantes foram organizados em uma tabela em nível de amostra e integrados ao metadata, gerando o arquivo `results/alpha_diversidade_por_amostra.tsv`, que pode ser utilizado em análises estatísticas subsequentes e em visualizações gráficas.[cite:563][cite:564]

### Agrupamento e visualização

Os valores de diversidade alfa por amostra foram visualizados por meio de boxplots agrupados pela coluna do metadata correspondente ao supersítio. Os rótulos dos gráficos foram padronizados em inglês para facilitar o reaproveitamento das figuras em manuscritos, apresentações e relatórios.[cite:725]

Foram gerados dois gráficos principais:

- diversidade de Shannon por supersítio.[cite:727]
- riqueza observada de ASVs por supersítio.[cite:725][cite:726]

### Exportação das figuras

As figuras de diversidade alfa foram exportadas nos formatos PNG e PDF com a função `ggplot2::ggsave`, amplamente utilizada para salvar objetos do `ggplot2` com definição explícita de nome do arquivo, dimensões e formato de saída.[cite:714][cite:733]

Os arquivos foram salvos no diretório:

```r
results/figures/
```

com os seguintes nomes:

```r
alpha_shannon_by_supersite.png
alpha_shannon_by_supersite.pdf
alpha_richness_by_supersite.png
alpha_richness_by_supersite.pdf
```

### Notas de reprodutibilidade

Esse fluxo depende da harmonização prévia dos identificadores de ASVs entre a tabela de abundância fúngica e a tabela expandida de taxonomia, bem como da correspondência exata entre os IDs das amostras na matriz de abundância e no arquivo `metadata.tsv`. A verificação `all(sample_names %in% rownames(meta_tab))` foi utilizada para confirmar que todas as amostras presentes na matriz de abundância estavam representadas no metadata antes do cálculo das métricas de diversidade.[cite:563][cite:574]

Os produtos gerados fornecem uma base reprodutível para comparar a diversidade fúngica intra-amostra entre supersítios e para análises inferenciais subsequentes, como testes não paramétricos ou correlações com variáveis ambientais.[cite:732][cite:726]

## Análise de diversidade beta por Bray-Curtis

A diversidade beta foi analisada em nível de amostra a partir da tabela de abundância de ASVs fúngicos com taxonomia (`results/asv_table_fungi_with_taxonomy.tsv`) e do arquivo de metadados (`metadata.tsv`). O objetivo foi quantificar diferenças na composição das comunidades fúngicas entre amostras e avaliar se os supersítios apresentam estrutura comunitária distinta.[cite:772][cite:778]

### Tabelas de entrada e organização dos dados

A tabela `results/asv_table_fungi_with_taxonomy.tsv` foi lida em R e separada em duas partes: uma seção taxonômica e uma matriz de abundância. A coluna `Confidence` foi utilizada como referência para identificar o ponto de separação entre taxonomia e contagens.[cite:772]

A matriz de abundância foi então transposta para o formato amostras x ASVs, que é a estrutura esperada para o cálculo de dissimilaridades ecológicas no pacote `vegan`.[cite:772][cite:761]

O arquivo `metadata.tsv` foi lido com os identificadores das amostras definidos como nomes de linha, e a compatibilidade entre os nomes das amostras da matriz de abundância e do metadata foi verificada antes das análises.[cite:563][cite:574]

### Matriz de dissimilaridade Bray-Curtis

A dissimilaridade entre amostras foi calculada com a função `vegdist(..., method = "bray")` do pacote `vegan`. A distância de Bray-Curtis utiliza a abundância das ASVs e é amplamente empregada em estudos de microbioma para comparar diferenças na estrutura das comunidades entre amostras.[cite:761][cite:782]

A matriz resultante foi utilizada como base tanto para a ordenação multivariada quanto para o teste estatístico entre grupos.[cite:777][cite:778]

### Ordenação por NMDS

A representação bidimensional da variação composicional foi obtida por NMDS (Non-metric Multidimensional Scaling) com a função `metaMDS`, utilizando duas dimensões (`k = 2`), múltiplas tentativas de convergência (`trymax = 100`) e sem transformação automática adicional dos dados (`autotransform = FALSE`). O algoritmo `metaMDS` busca uma solução estável por múltiplos inícios aleatórios, o que melhora a robustez da ordenação.[cite:773][cite:781]

As coordenadas das amostras foram extraídas com `scores(..., display = "sites")` e integradas ao metadata para permitir a visualização dos pontos coloridos por supersítio em `ggplot2`.[cite:775][cite:765]

### Teste entre supersítios por PERMANOVA

A significância das diferenças na composição fúngica entre supersítios foi avaliada por PERMANOVA com a função `adonis2`, utilizando a matriz de Bray-Curtis como resposta e o fator supersítio como variável explicativa. O teste foi executado com 999 permutações.[cite:777][cite:767]

A saída do `adonis2` fornece, para o fator testado, a soma de quadrados, a estatística pseudo-\(F\), a proporção de variação explicada (`R2`) e o valor de significância permutacional (`Pr(>F)`). Esses parâmetros permitem avaliar se a composição fúngica difere significativamente entre supersítios e qual fração da variação total é explicada por esse agrupamento.[cite:777][cite:651]

### Visualização e exportação das figuras

A ordenação NMDS foi representada graficamente em `ggplot2`, com pontos correspondentes às amostras e coloração de acordo com o supersítio. Os rótulos foram padronizados em inglês para facilitar o uso em manuscritos e apresentações.[cite:765]

A figura foi exportada em formatos PNG e PDF com `ggplot2::ggsave`, definindo-se explicitamente o nome do arquivo, as dimensões e, no caso do PNG, a resolução em dpi.[cite:708][cite:714]

Os arquivos foram destinados ao diretório:

```r
results/figures/
```

com os nomes:

```r
nmds_bray_by_supersite.png
nmds_bray_by_supersite.pdf
```

### Reprodutibilidade

Esse fluxo depende da correspondência exata entre os identificadores das amostras na tabela de abundância e no arquivo `metadata.tsv`, bem como da correta especificação da coluna do metadata utilizada para representar o supersítio. A reordenação explícita do metadata com base em `sample_names` foi utilizada para garantir consistência entre a matriz comunidade e os fatores amostrais empregados na PERMANOVA e na visualização.[cite:563][cite:574]

Os produtos gerados por essa etapa constituem a base para interpretar padrões de beta diversidade entre supersítios e podem ser complementados, em etapas posteriores, por análises com Jaccard, testes de dispersão multivariada e figuras com elipses ou envoltórias por grupo.[cite:647][cite:782]

## Avaliação da dispersão multivariada por Bray-Curtis

Como complemento à PERMANOVA baseada em Bray-Curtis, foi realizada uma análise de homogeneidade de dispersões multivariadas entre supersítios por meio da função `betadisper` do pacote `vegan`. Essa abordagem corresponde ao procedimento PERMDISP, utilizado para avaliar se os grupos diferem em sua variabilidade interna no espaço multivariado.[cite:788][cite:789]

### Estrutura de entrada

A análise utilizou como entrada a matriz de dissimilaridade Bray-Curtis previamente calculada entre amostras (`dist_bray`) e o fator de agrupamento correspondente ao supersítio, extraído do metadata reordenado (`meta_tab2`). A correspondência entre amostras e grupos foi mantida pela mesma estrutura empregada na NMDS e na PERMANOVA anteriores.[cite:788][cite:789]

### Estimativa da dispersão por grupo

A função `betadisper(dist_bray, group = group_vec)` foi empregada para calcular a distância de cada amostra ao centro do seu respectivo grupo no espaço derivado da matriz de dissimilaridade. No `vegan`, esse procedimento é tratado como um análogo multivariado do teste de Levene para homogeneidade de variâncias.[cite:788]

As distâncias ao centróide ou à mediana do grupo podem ser usadas para quantificar dispersão interna; neste fluxo, a saída foi utilizada para resumir a distância média por supersítio e para derivar os valores individuais empregados na visualização por boxplot.[cite:788][cite:789]

### Teste permutacional

A significância das diferenças entre dispersões foi testada com `permutest(bd_bray, permutations = 999)`, utilizando 999 permutações livres. Esse procedimento gera uma distribuição nula para a estatística \(F\) sob a hipótese de ausência de diferença na dispersão entre grupos.[cite:789][cite:799]

A interpretação do teste considera que valores significativos indicam heterogeneidade de dispersão entre os supersítios, o que é relevante para contextualizar a leitura dos resultados de PERMANOVA obtidos a partir da mesma matriz de dissimilaridade.[cite:801][cite:807]

### Visualização gráfica

As distâncias individuais ao centróide foram organizadas em um `data.frame` com uma coluna para o supersítio e outra para a distância observada (`bd_bray$distances`). Esses valores foram representados graficamente com `ggplot2` em boxplots por grupo, permitindo comparar a variação interna das amostras em cada supersítio.[cite:788]

Os eixos foram rotulados em inglês, com o objetivo de facilitar o reaproveitamento da figura em manuscritos, relatórios e apresentações. A figura foi salva em PNG e PDF com `ggsave`, com dimensões explícitas para padronização visual dos resultados.[cite:714][cite:733]

Os arquivos de saída previstos para essa etapa foram:

```r
results/figures/betadisper_bray_by_supersite.png
results/figures/betadisper_bray_by_supersite.pdf
```

### Papel analítico no estudo

A análise de PERMDISP foi empregada para verificar se diferenças detectadas por PERMANOVA entre supersítios poderiam estar associadas, ao menos em parte, a diferenças na dispersão interna dos grupos. Assim, essa etapa funciona como verificação complementar na interpretação da diversidade beta baseada em Bray-Curtis.[cite:789][cite:807]

## Análise de beta-diversidade (Supersitio)

A matriz de abundância de ASVs fúngicas por amostra (`ASV_abundance_fungi.tsv`) foi convertida em presença/ausência, atribuindo valor 1 para abundâncias maiores que zero e 0 para ausências. Essa binarização foi feita linha a linha (por amostra), de forma que a análise ficasse baseada apenas na composição de ASVs, sem levar em conta diferenças de abundância relativa. [web:833][file:591]

A partir da matriz binária por amostra, foi calculada a dissimilaridade entre amostras usando o índice de Jaccard no pacote **vegan**, via `vegdist(..., method = "jaccard", binary = TRUE)`. No **vegan**, o argumento `binary = TRUE` é a forma recomendada para tratar explicitamente a matriz como presença/ausência ao usar o índice de Jaccard. [web:761][web:833]

A estrutura de dissimilaridade da comunidade fúngica entre amostras foi então resumida por meio de escalonamento multidimensional não métrico (NMDS), utilizando a função `metaMDS` com duas dimensões (`k = 2`) e distância de Jaccard binária como entrada. O ajuste da ordenação foi avaliado pelo valor de *stress* retornado pelo `metaMDS`, que reflete o quanto a configuração bidimensional preserva as relações de dissimilaridade entre as amostras. [web:763][web:827]

Para testar diferenças de composição da comunidade entre supersítios, foi aplicada uma PERMANOVA com a função `adonis2` (pacote **vegan**), usando a matriz de distâncias de Jaccard binária e o fator `Supersitio` definido no arquivo `metadata.tsv`. Foram utilizadas 999 permutações, e a estatística F, o R² e o valor de p foram usados para interpretar a significância e o tamanho de efeito do fator. Esse procedimento segue recomendações correntes para análises multivariadas de dados de metabarcoding em ecologia de comunidades. [web:767][file:591]

A figura de ordenação foi construída a partir dos escores do NMDS, com cada ponto representando uma amostra e as cores indicando os diferentes níveis de `Supersitio`. Elipses de dispersão foram adicionadas para cada supersítio, de forma a ilustrar a variabilidade interna dos grupos no espaço de ordenação. A legenda da figura inclui o valor de *stress* do NMDS e o resumo da PERMANOVA (F, R² e p para `Supersitio`). [web:827][file:591]

## Beta-diversity analysis by supersite

The fungal ASV abundance table was imported from `results/qiime2/feature-table_fungi_export/ASV_abundance_fungi.tsv`, with ASVs arranged in rows and samples in columns. The table was transposed so that samples represented rows and ASVs represented columns, which is the format required for downstream community dissimilarity analyses in `vegan`. [web:772][web:874]

To focus the analysis on community composition rather than abundance, the abundance matrix was converted to a binary presence/absence matrix. For each ASV in each sample, values greater than zero were recoded as 1 and zero values were kept as 0. This approach is appropriate when the objective is to evaluate ASV turnover and compositional similarity among samples independently of relative abundance. [web:870][web:819]

Pairwise dissimilarity among samples was calculated using the Jaccard index with `vegdist(..., method = "jaccard", binary = TRUE)` in the `vegan` package. In this implementation, `binary = TRUE` ensures that the distance is computed from presence/absence data rather than abundance values. [web:819][web:791]

Ordination was performed using non-metric multidimensional scaling (NMDS) with the function `metaMDS`, using the binary Jaccard distance matrix as input and a two-dimensional solution (`k = 2`). The fit of the ordination was assessed by the stress value returned by the model, which quantifies how well the two-dimensional representation preserves the rank-order structure of the original dissimilarity matrix. [web:864][web:772]

Differences in fungal community composition among supersites were tested using permutational multivariate analysis of variance (PERMANOVA) with the function `adonis2`, using 999 permutations and `Supersitio` as the explanatory factor from `metadata.tsv`. The PERMANOVA output was interpreted using the pseudo-F statistic, the proportion of explained variation (R²), and the permutation-based p-value. [web:864][web:874]

The NMDS figure was generated with `ggplot2`, with points representing samples and colors indicating the levels of `Supersitio`. Group ellipses were added to summarize within-group dispersion in ordination space. For reproducibility and downstream use in manuscripts and presentations, all figures were exported in both PNG and PDF formats using `ggsave()`. [web:708][web:714]

## Análise de beta-diversidade por supersítio

A tabela de abundância de ASVs fúngicas foi importada a partir do arquivo `results/qiime2/feature-table_fungi_export/ASV_abundance_fungi.tsv`, no qual as ASVs estavam organizadas em linhas e as amostras em colunas. Para as análises de ecologia de comunidades, a matriz foi transposta, de modo que as amostras passassem a ocupar as linhas e as ASVs as colunas. [web:772][web:911]

Com o objetivo de avaliar a composição da comunidade independentemente da abundância, a matriz de abundância foi convertida em uma matriz binária de presença/ausência. Valores maiores que zero foram recodificados como 1, enquanto valores iguais a zero foram mantidos como 0. Essa abordagem permite que a análise reflita exclusivamente o compartilhamento e a substituição de ASVs entre amostras. [web:877][web:911]

A dissimilaridade entre amostras foi calculada com o índice de Jaccard utilizando a função `vegdist(..., method = "jaccard", binary = TRUE)` do pacote `vegan`. Nessa implementação, o argumento `binary = TRUE` garante que a distância seja calculada a partir de dados de presença/ausência. [web:772][web:911]

A ordenação das amostras foi realizada por escalonamento multidimensional não métrico (NMDS), com a função `metaMDS`, utilizando uma solução bidimensional (`k = 2`) baseada na matriz de dissimilaridade de Jaccard binária. A qualidade do ajuste da ordenação foi avaliada pelo valor de *stress*, que expressa o quanto a representação bidimensional preserva a estrutura de distâncias da matriz original. [web:772][web:911]

As diferenças na composição da comunidade fúngica entre supersítios foram testadas por PERMANOVA com a função `adonis2`, utilizando `Supersitio` como fator explicativo e 999 permutações. A interpretação do teste foi baseada na estatística pseudo-F, na proporção de variação explicada (R²) e no valor de p obtido por permutação. [web:913][web:914]

Como resultados significativos de PERMANOVA podem ser influenciados por diferenças na dispersão interna dos grupos, a homogeneidade das dispersões multivariadas foi avaliada com a função `betadisper`, usando a mesma matriz de Jaccard binária e o fator `Supersitio`. Esse procedimento constitui um análogo multivariado do teste de Levene e estima a distância média das amostras ao centróide do grupo. A significância foi avaliada tanto por ANOVA quanto por teste de permutação com `permutest(..., permutations = 999)`. [web:877][web:800][web:788]

Foram geradas duas figuras para esta análise: um gráfico de ordenação NMDS e um boxplot das distâncias ao centróide derivadas do `betadisper`. Para fins de reprodutibilidade e uso em manuscritos, todas as figuras foram exportadas em formatos PNG e PDF com a função `ggsave()`. [web:910][web:915]