# Nome do projeto
<!-- Aqui o aluno coloca um título simples e direto.
     Ex.: "QC de dados de RNA-Seq" ou "Análise de FASTQ de orquídeas". -->

Breve descrição do projeto em 1 ou 2 frases.
<!-- Dizer em português claro o que o projeto faz.
     Ex.: "Este projeto faz o controle de qualidade de arquivos FASTQ usando FastQC e MultiQC." -->

## Objetivo
<!-- Seção para explicar o PORQUÊ do projeto existe.
     Falar em linguagem simples, sem muito jargão. -->

Exemplo:
- Organizar um projeto de bioinformática em pastas.
- Rodar ferramentas de controle de qualidade em arquivos FASTQ.
- Gerar relatórios para avaliar a qualidade dos dados de sequenciamento.

<!-- O aluno pode adaptar esses bullets para o objetivo real do projeto dele. -->

## Estrutura do projeto
<!-- Aqui o aluno mostra COMO os arquivos estão organizados.
     Isso ajuda qualquer pessoa a se achar dentro do projeto. -->

```text
data/raw/         # arquivos de entrada (ex.: FASTQ)
data/metadata/    # informações das amostras (ex.: samples.csv)
scripts/          # scripts do projeto (bash, python, R)
results/          # resultados gerados (relatórios, tabelas, figuras)
README.md         # este arquivo de documentação
```

<!-- O aluno pode acrescentar ou remover pastas conforme o projeto,
     mas deve manter os comentários explicando o uso de cada uma. -->

## Programas necessários
<!-- Lista de tudo que precisa estar instalado ANTES de rodar o projeto.
     Serve como checklist para outra pessoa reproduzir o ambiente. -->

Exemplo:
- VS Code (para editar códigos e abrir o terminal)
- Conda ou Mamba (para gerenciar ambientes)
- Git (para controle de versão, opcional na primeira vez)
- FastQC
- MultiQC

<!-- O aluno pode trocar ou adicionar programas (por ex.: Trimmomatic, STAR, R, etc.). -->

## Como criar o ambiente
<!-- Aqui entram os comandos que a pessoa deve rodar no terminal
     para preparar o ambiente (conda, pacotes, etc.). -->

```bash
# criar ambiente (o nome pode ser adaptado)
conda create --name meu-projeto python=3.11 -y

# ativar ambiente
conda activate meu-projeto

# instalar ferramentas necessárias
conda install fastqc multiqc -y
```

<!-- Orientar o aluno a testar se deu certo com:
     fastqc --help
     multiqc --help
-->

## Como rodar
<!-- Seção fundamental: passo a passo de execução.
     Pense em alguém que NUNCA viu o projeto e só vai copiar e colar. -->

Passo 1: ativar o ambiente

```bash
conda activate meu-projeto
```

Passo 2: rodar o(s) script(s) principal(is)

```bash
./scripts/run_fastqc.sh
./scripts/run_multiqc.sh
```

<!-- O aluno deve adaptar esses comandos para o que ele realmente usa.
     Ex.: ./scripts/run_pipeline.sh, python scripts/analise.py, etc. -->

## Dados de entrada
<!-- Explicar ONDE colocar os arquivos que serão analisados
     e qual o formato esperado (nome, extensão, tipo de dado). -->

Exemplo:
- Coloque os arquivos FASTQ na pasta `data/raw/`.
- Os arquivos podem terminar em `.fastq` ou `.fastq.gz`.
- Exemplo de nomes: `amostra1.fastq.gz`, `amostra2.fastq.gz`.

<!-- Se houver tabela de metadados, citar aqui:
     "O arquivo data/metadata/samples.csv contém as informações das amostras." -->

## Resultados
<!-- Explicar ONDE os resultados aparecem depois que o projeto roda
     e QUAL o tipo de saída (HTML, TXT, figuras, etc.). -->

Exemplo:
- Relatórios individuais do FastQC: pasta `results/fastqc/`.
- Relatório consolidado do MultiQC: pasta `results/multiqc/`
  (arquivo `multiqc_report.html`).

<!-- O aluno pode acrescentar outros resultados, como:
     - "Tabela de contagens em results/counts/"
     - "Figuras em results/figures/" -->

## Autor
<!-- Espaço para o aluno se identificar. -->

Nome do aluno ou grupo.  
Disciplina / Curso / Instituição (opcional).

## Observações
<!-- Campo livre para anotar qualquer detalhe importante:
     limitações, pressupostos, ideias futuras, etc. -->

Exemplos:
- Este projeto foi desenvolvido como atividade prática de bioinformática.
- Pode ser usado como modelo para outros projetos mudando apenas os arquivos de entrada.
