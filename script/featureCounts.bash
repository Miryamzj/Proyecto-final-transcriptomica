set -euo pipefail
# ==========================================
# featureCounts
# Proyecto: Zombie ant RNA-seq
# ==========================================
PROJECT=$HOME/work_directory/zombie_ant_rnaseq
FEATURECOUNTS="/export/apps/bioconda20240614/envs/subread/bin/featureCounts"
GTF="$PROJECT/reference/GCF_003227725.1_Cflo_v7.5_genomic.gtf"
THREADS=10
BAMDIR="$PROJECT/results/star"
OUTDIR="$PROJECT/results/featureCounts"
mkdir -p "$OUTDIR"
echo "Iniciando featureCounts..."
date
/usr/bin/time -v "$FEATURECOUNTS" \-T ${THREADS} \ # Número de hilos-a "$GTF" \ # Archivo de anotación GTF-o "$OUTDIR/gene_counts.txt" \ # Archivo de salida con conteos-t exon \ # Tipo de característica a contar-g gene_id \ # Atributo de identificación génica--largestOverlap \ # Asignar a la característica con mayor solapamiento-p \ # Contar fragmentos en lugar de lecturas-P \ # Considerar solo pares correctamente emparejados-B \ # Requerir que ambos extremos estén alineados
${BAMDIR}/*.bam \ # Archivos BAM de entrada
> "$OUTDIR/featureCounts.stdout.log" \ # Log de salida estándar
2> "$OUTDIR/featureCounts.time.log"  # Log de errores y tiempo
echo "Finalizado."
date