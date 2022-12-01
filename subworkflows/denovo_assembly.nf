import groovy.json.JsonSlurper
import nextflow.util.BlankSeparatedList;


map_sample_ids_cls = {it -> 
/* Harmonize tuples
output:
    tuple val(sample_id), path('*.gff')
When there are multiple paths, will emit:
    [sample_id, [path, path ..]]
when there's a single path, this:
    [sample_id, path]
This closure makes both cases:
    [[sample_id, path][sample_id, path]].
*/
    if (it[1].getClass() != java.util.ArrayList){
        // If only one path, `it` will be [sample_id, path]
        return [it]
    }
    l = [];
    for (x in it[1]){
        l.add(tuple(it[0], x))
    }
    return l
}


process dump_clusters {
    label "isoforms"
    input:
        tuple val(sample_id), path(root_cluster), path(sorted_reads_dir)
    output:
        tuple val(sample_id), path("final_clusters"), emit: final_clusters_dir
        tuple val(sample_id), path("final_clusters/cluster_fastq/*.fq"), emit: final_clusters
    shell:
    """ isONclust2 dump -v -i $sorted_reads_dir/sorted_reads_idx.cer -o final_clusters $root_cluster; sync """

}

process build_backbones {
    /*
    Build coding
    */
    label "isoforms"
	input:
		tuple val(sample_id), path(cluster_fq)
	output:
		tuple val(sample_id), path("*final_polished_cds.fa"), emit: polished_cds, optional: true
	script:
	def cluster_fq_bl = new BlankSeparatedList(cluster_fq)
	"""
    # Get one of the cluster ids to give the output a unique name
    UNID=\$(echo ${cluster_fq_bl[1]} | grep -o -E '[0-9]+')

	for cluster in $cluster_fq_bl
        do
            clfq=`basename \$cluster`
            cln=\${clfq%.*}

            echo Building backbone for cluster: \$cln
            echo "\tSampling input reads for backbone construction."
            sample=\${cln}_sample.fq
            seqkit head --quiet -n 100 \$cluster > \$sample
            seqkit sample --quiet -n 500 -2 -s 100 \$cluster >> \$sample
            echo "\tConstructing spoa consensus."
            spoa_cons=\${cln}_spoa.fa
            spoa -m 5 -n -4 -g -8 -e -6 -q -10 -c -15 -l 1 -r 0 \$sample > \$spoa_cons

            echo "\tPolishing the consensus using racon."

            # polish 1
            samgz=\${cln}_aln.sam.gz
            racon_cons=\${cln}_racon.fa
            tmpcons=\${cln}_tmcons.fa

            minimap2 -t ${params.threads} -ax splice ${params.minimap2_opts} \$spoa_cons \$sample | gzip - > \$samgz
            exitcode=0
            racon -t ${params.threads} --no-trimming -u -w 2000 \$sample \$samgz \$spoa_cons > \$racon_cons || exitcode=1
            if [[ \$exitcode -eq 0 ]];
                then
                    # Rename consensus sequence name with cluster id
                    cat \$racon_cons | seqkit replace -p ".*" -r cluster_\${cln} > \$tmpcons
                    mv \$tmpcons \$racon_cons
            else
                echo "Polishing failed for \${cln}"
                continue
            fi

            # polish 2
            minimap2 -t ${params.threads} -ax splice ${params.minimap2_opts} \$racon_cons \$sample  | gzip - > \$samgz
            exitcode=0
            racon -t ${params.threads}  -u --no-trimming \$sample \$samgz \$racon_cons > \$tmpcons || exitcode=1
            if [[ \$exitcode -eq 0 ]];
                then
                    echo "success polish 2"
                    mv \$tmpcons \$racon_cons
            else
                echo "Polishing failed for \${cln}"
                continue
            fi

             # polish 3
            minimap2 -t ${params.threads} -ax splice ${params.minimap2_opts} \$racon_cons \$sample | gzip - > \$samgz
            exitcode=0
            racon -t ${params.threads}  -u \$sample \$samgz \$racon_cons > \$tmpcons || exitcode=1
            if [[ \$exitcode -eq 0 ]];
                then
                    cat \$tmpcons >> ${sample_id}_\${UNID}_final_polished_cds.fa
                    echo "polishing 3 success"
            else
                echo "Polishing failed for \${cln}"
            fi

        done
    echo "Finished backbones"
	"""
}


process merge_cds {
    label "isoforms"
    input:
        tuple val(sample_id), path(cds)
    output:
        tuple val(sample_id), path("${sample_id}_cds.fa"), emit: final_polished_cds
    script:
    """
    for FILE in *final_polished_cds.fa
    do
        cat \$FILE  >> "${sample_id}_cds.fa"
    done
    """
}

process cds_align {
    label "isoforms"
    input:
        tuple val(sample_id), path(polished_cds), path(sorted_reads_dir)
    output:
        tuple val(sample_id), path("${sample_id}_reads_aln_sorted.bam"), emit: bam
        tuple val(sample_id), path("${sample_id}_read_aln_stats.tsv"), emit: stats
    script:
    """
    minimap2 -t ${params.threads} \
        -ax splice ${params.minimap2_opts} $polished_cds $sorted_reads_dir/sorted_reads.fastq |\
        samtools view -b - |\
    samtools sort -o "${sample_id}_reads_aln_sorted.bam";
    samtools index "${sample_id}_reads_aln_sorted.bam";
    ((seqkit bam -s -j ${params.threads}  "${sample_id}_reads_aln_sorted.bam" 2>&1)  | tee ${sample_id}_read_aln_stats.tsv ) || true
	"""
}


process make_batches {
    /*
    Take a fasta file and creates batches for isONclust2 to work with

    */
    label "isoforms"
    input:
        tuple val(sample_id), path(fastq)
    output:
        tuple val(sample_id), path('sorted/batches'), emit: sorted_batches
        tuple val(sample_id), path('sorted'), emit: sorted_reads_dir
    
    script:

    maxcpus = Runtime.runtime.availableProcessors()

    minimum_batch_size = 2000
    """
    b=0
    if [ ${params.isOnClust2_batch_size} -lt \$b ];
    then
        nr_bases=\$(seqkit stats -T $fastq|cut -f 5| sed '2q;d')
        let batch_size=\$nr_bases/1000/$maxcpus
        if [ \$batch_size -lt $minimum_batch_size ];
            then
                batch_size=$minimum_batch_size
        fi
    else
        batch_size=${params.isOnClust2_batch_size}
    fi

    echo "Batch size:\$batch_size";
    echo "Num bases: \$nr_bases";

    mkdir -p sorted; isONclust2 sort $params.isOnClust2_sort_options -v -o sorted $fastq;
    """

}

process clustering() {
    label "isoforms"

    input:
        tuple val(sample_id), path(sorted_batches)
    output:
        tuple val(sample_id), path('isONcluster_ROOT.cer'), emit: root_cluster
    script:
    """
    run_isonclust2.py $sorted_batches
    """
}

process cluster_quality() {
    // Run Kristoffer Sahlin's QC code.
    // For now just write out the PDF and CSV results. Move these to the report at some point
    label "isoforms"

    input:
        path reference
        tuple val(sample_id), path(reads_fl), path(final_clusters_dir)
    output:
        tuple val(sample_id),
            path("${sample_id}_cluster_qc"), emit: cluster_qc_dir
        tuple val(sample_id),
            path("${sample_id}_cluster_qc_raw"), emit: cluster_qc_raw
    script:
    def qc_dir = "${sample_id}_cluster_qc"
    def qc_dir_raw = "${sample_id}_cluster_qc_raw" // To generate plots in report
    def bam = "${qc_dir}/ref_aln.bam"

    """
    mkdir $qc_dir
    mkdir $qc_dir_raw
    minimap2 -ax splice -t 2 $reference $reads_fl |\
        samtools view -q 2 -F 2304 -b - |\
         samtools sort - -o $bam;
    samtools index $bam;
    compute_cluster_quality.py --sizes $final_clusters_dir/clusters_info.tsv \
        --outfile ${qc_dir}/cluster_quality.csv --ont --clusters $final_clusters_dir/clusters.tsv \
        --classes $bam --report ${qc_dir}/cluster_quality.pdf --raw_data_out $qc_dir_raw
    """

}

workflow denovo_assembly {
    take:
       fastq_reads_fl
       reference
    main:
        make_batches(fastq_reads_fl)

        clustering(make_batches.out.sorted_batches)

        dump_clusters(
            clustering.out.root_cluster
            .join(make_batches.out.sorted_reads_dir))

        build_backbones(
            dump_clusters.out.final_clusters
            .flatMap(map_sample_ids_cls)
            .groupTuple(size: 10, remainder: true))

        merge_cds(
            build_backbones.out.polished_cds
            .flatMap(map_sample_ids_cls)
            .groupTuple())

        cds_align(
            merge_cds.out.final_polished_cds
            .join(make_batches.out.sorted_reads_dir))

        if (!reference.name.startsWith('OPTIONAL_FILE')){
            cluster_quality(reference, fastq_reads_fl
            .join(dump_clusters.out.final_clusters_dir))

            cluster_quality.out.cluster_qc_dir.set { opt_qual_ch }

            cluster_quality.out.cluster_qc_raw.set { opt_qual_raw_ch }

        } else{
            Channel.empty().set { opt_qual_ch }
            Channel.empty().set { opt_qual_raw_ch }
        }

    emit:
       bam = cds_align.out.bam
       cds = merge_cds.out.final_polished_cds
       stats = cds_align.out.stats
       opt_qual_ch
       opt_qual_raw_ch
}

