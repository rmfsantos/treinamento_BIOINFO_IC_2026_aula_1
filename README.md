# Projeto de QC com FastQC e MultiQC

Este projeto foi criado para analisar a qualidade de arquivos FASTQ usando as ferramentas **FastQC** e **MultiQC**.

## Objetivo

Neste projeto, o aluno vai:

- Organizar um projeto de bioinformática em pastas.
- Rodar o FastQC para avaliar a qualidade dos arquivos FASTQ.
- Rodar o MultiQC para juntar todos os relatórios em um único arquivo.
- Aprender a documentar o projeto com um README.

## Estrutura do projeto

```text
data/raw/          # arquivos FASTQ de entrada
data/metadata/     # informações das amostras
scripts/           # scripts para rodar as análises
results/fastqc/    # resultados do FastQC
results/multiqc/   # resultado do MultiQC
README.md          # documentação do projeto
```

## Programas necessários

Antes de começar, é preciso ter instalado:

- Crhome
- Libre office
- VS Code
- Mamba
- Git (opcional, mas recomendado)

## Como criar o ambiente

Se quiser criar o ambiente manualmente, use:

```bash
conda create --name orquideas-qc python=3.11 -y
conda activate orquideas-qc
conda install fastqc multiqc -y #instala as ferramentas necessarias
```

## Dados de entrada

Os arquivos de entrada devem ser colocados em:

```text
data/raw/
```

Exemplo de arquivos:

```text
data/raw/amostra1.fastq.gz
data/raw/amostra2.fastq.gz
data/raw/amostra3.fastq.gz
```

Se o projeto tiver metadados, eles podem ser colocados em:

```text
data/metadata/samples.csv
```

## Script do FastQC

Arquivo: `scripts/run_fastqc.sh`

```bash
#!/usr/bin/env bash
# Define o interpretador do script usando o Bash localizado no ambiente do sistema

set -euo pipefail
# -e: interrompe o script se qualquer comando falhar
# -u: gera erro se uma variável não definida for usada
# -o pipefail: falha se qualquer comando em um pipeline falhar

RAW_DIR="data/raw"
# Diretório onde estão os arquivos FASTQ de entrada

OUT_DIR="results/fastqc"
# Diretório onde os resultados do FastQC serão armazenados

echo "Criando diretório de saída: ${OUT_DIR}..."
# Mensagem informando a criação do diretório de saída

mkdir -p "${OUT_DIR}"
# Cria o diretório de saída; -p evita erro se já existir

echo "Rodando FastQC em todos os arquivos FASTQ de ${RAW_DIR}..."
# Mensagem indicando início da análise

fastqc "${RAW_DIR}"/*.fastq* -o "${OUT_DIR}"
# Executa o FastQC em todos os arquivos .fastq, .fq, .fastq.gz etc.
# dentro do diretório RAW_DIR, salvando os resultados em OUT_DIR

echo "Análise FastQC concluída."
# Mensagem indicando término da execução

echo "Relatórios salvos em: ${OUT_DIR}"
# Informa onde os relatórios foram salvos
```

Esse padrão de script segue o uso típico do FastQC em lote sobre múltiplos arquivos FASTQ.

## Script do MultiQC

Arquivo: `scripts/run_multiqc.sh`

```bash
#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="orquideas-qc"
IN_DIR="results/fastqc"
OUT_DIR="results/multiqc"

echo "Criando diretório de saída do MultiQC: ${OUT_DIR}..."
mkdir -p "${OUT_DIR}"

echo "Rodando MultiQC em ${IN_DIR}..."
multiqc "${IN_DIR}" -o "${OUT_DIR}"

echo "Análise MultiQC concluída."
echo "Relatório principal: ${OUT_DIR}/multiqc_report.html"
```

O MultiQC é usado justamente para agregar vários relatórios de ferramentas como FastQC em um único HTML.

## Como rodar o projeto

### 1. Ativar o ambiente

```bash
conda activate orquideas-qc
```

### 2. Dar permissão aos scripts

```bash
chmod +x scripts/run_fastqc.sh
chmod +x scripts/run_multiqc.sh
```

### 3. Rodar o FastQC

```bash
./scripts/run_fastqc.sh
```

### 4. Rodar o MultiQC

```bash
./scripts/run_multiqc.sh
```

Esse fluxo corresponde ao procedimento básico ensinado em tutoriais introdutórios de QC de reads com FastQC e MultiQC.[web:22][web:64][web:11]  

## Resultados

Depois da execução:

- Os relatórios do FastQC ficarão em:
  - `results/fastqc/`
- O relatório final do MultiQC ficará em:
  - `results/multiqc/multiqc_report.html`

O relatório MultiQC permite comparar várias amostras de uma vez, o que facilita a interpretação do QC.[web:11][web:80][web:81]  

## Autor

Preencha com:

- Nome do aluno ou grupo
- Disciplina
- Instituição

Exemplo:

- Aluno: Maria Silva
- Disciplina: Introdução à Bioinformática
- Instituição: UEFS

## Observações

Este projeto pode ser reutilizado em outros conjuntos de dados, bastando trocar os arquivos FASTQ na pasta `data/raw/` e ajustar os metadados, se necessário.
