OUTDIR=$PROJECT/results/star
LOGDIR=$PROJECT/results/logs

mkdir -p $OUTDIR
mkdir -p $LOGDIR

for R1 in ${FASTQ_DIR}/*_1.fastq
do

    SAMPLE=$(basename $R1 _1.fastq)
    R2=${FASTQ_DIR}/${SAMPLE}_2.fastq

    echo "===================================="
    echo "Procesando: $SAMPLE"
    echo "===================================="

    BAM=${OUTDIR}/${SAMPLE}.bam

    if [ -f "$BAM" ]; then
        echo "$SAMPLE ya existe. Saltando..."
        continue
    fi

    STAR \
        --runThreadN ${THREADS} \ # numero de hilos a usar
        --genomeDir ${STAR_INDEX} \ # directorio del indice de STAR
        --readFilesIn ${R1} ${R2} \ # archivos de entrada
        --outSAMtype BAM SortedByCoordinate \ # formato de salida y ordenamiento BAM
        --quantMode GeneCounts \ # habilitar conteo de genes
        --outFileNamePrefix ${OUTDIR}/${SAMPLE}_ \ # prefijo para los archivos de salida
        > ${LOGDIR}/${SAMPLE}.log 2>&1 # redirigir stdout y stderr al archivo de log

    mv ${OUTDIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam $BAM

    samtools index $BAM

    samtools flagstat $BAM \
        > ${OUTDIR}/${SAMPLE}.flagstat.txt

done
