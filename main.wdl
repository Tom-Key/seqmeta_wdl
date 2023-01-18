version 1.0


workflow wf_seqmeta{
    input {
        String sample_id
        Array[File] sample_fq1 
        Array[File] sample_fq2 
        
        String processdir
        String docker_img = "registry.cn-shenzhen.aliyuncs.com/oseq-rd/seqmeta:0.6"
    }


    call cat {
        input:
            MGISeqData_fq1 = sample_fq1, 
            MGISeqData_fq2 = sample_fq2,
            sample_id = sample_id ,
            processdir = processdir,
            docker_img = docker_img
    }
    call fastp_pair_end{
        input:
            in = cat.out_fq,
            sample_id = sample_id ,
            processdir = processdir,
            docker_img = docker_img
    }

    call bowite2_align{
        input:
            fastq1 = fastp_pair_end.out1,
            fastq2 = fastp_pair_end.out2,
            sample_id = sample_id ,
            processdir = processdir,
            docker_img = docker_img
    }

    call hpv{
        input:
            rmhost_fastq_se = bowite2_align.rmhost_fastq_se,
            sample_id = sample_id ,
            processdir = processdir,
            docker_img = docker_img
    }

    call humann2{
        input:
            rmhost_fastq = bowite2_align.rmhost_fastq,
            sample_id = sample_id ,
            processdir = processdir,
            docker_img = docker_img
    }

    output{
        Array[File] bge_profile =  [
            humann2.humann2_result[0], humann2.humann2_result[1], humann2.humann2_result[2], humann2.metaphlan_bugs, humann2.marker_ab,
            hpv.hpv_out, hpv.seqkit_result, hpv.summary_L1,
            fastp_pair_end.log, bowite2_align.log, humann2.log,

        ]
    }
}

task cat{
    input{
        String sample_id
        Array[File] MGISeqData_fq1
        Array[File] MGISeqData_fq2
        String processdir        
        # Resource
        Int cpu = 2
        String memory = "4G"
        String disks = "local-disk 40 cloud_ssd"
        
        # docker image
        String docker_img = docker_img
    }
    
    String directory =  processdir+"1.assay/01.trimming/"
    String outname1 = directory+sample_id+".1.fq.gz"
    String outname2 = directory+sample_id+".2.fq.gz"
    
    command <<<
        mkdir -p ~{directory}
        cat ~{sep=" " MGISeqData_fq1}  > ~{outname1}
        cat ~{sep=" " MGISeqData_fq2}  > ~{outname2}
    >>>
    
    runtime {
        cpu: cpu
        memory: memory
        disks: disks
        docker: docker_img
    }
    output {
        Pair[File,File] out_fq = (outname1, outname2)
    }
}

task fastp_pair_end {

    input {
        
        # I/O options
        Pair [File,File] in

        Boolean? phred64 = false 
        Boolean? fix_mgi_id = false

        String? adapter_sequence
        String? adapter_sequence_r2

        Int? reads_to_process # specify how many reads/pairs to be processed. Default 0 means process all reads.
        
        String report_title = "\'fastp report\'"

        # excute env
        Int cpu = 4
        String memory = "8G"
        String disks = "local-disk 40 cloud_ssd"
        
        String processdir
        String sample_id
        String docker_img

    }


    String directory =  processdir+"1.assay/01.trimming/"
    String out1_name = directory+sample_id+".raw.1.fq.gz"
    String out2_name = directory+sample_id+".raw.2.fq.gz"
    # reporting options
    String json = directory+sample_id+".fastp.json"
    String html = directory+sample_id+".fastp.html"
    String log_directory =  processdir+"1.assay/logs/01.trimming/"
    String log = log_directory+sample_id+".fastp.log"
    

    command <<<
        source /root/miniconda3/bin/activate
        conda activate seqmeta
        cpu_cores=$(nproc)
        mkdir -p ~{directory}
        mkdir -p ~{log_directory}
        # basic command
        fastp \
        --in1 ~{in.left} \
        --in2 ~{in.right} \
        --out1 ~{out1_name} \
        --out2 ~{out2_name} \
        -w ~{cpu} \
        --compression 6 \
        --cut_front --cut_front_window_size 4 \
        --cut_front_mean_quality 20 \
        --cut_right \
        --cut_right_window_size 4 \
        --cut_right_mean_quality 20 \
        --n_base_limit 5 \
        --length_required 30 \
        --adapter_sequence=AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
        --adapter_sequence_r2=AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
        --json ~{json} \
        --html ~{html} \
        --report_title ~{report_title} 2>~{log}
    >>>

    runtime {
        cpu: cpu
        memory: memory
        disks: disks
        docker: docker_img
    }

    output {
        File out1 = out1_name
        File out2 = out2_name
        File json_report = json
        File html_report = html
        File log = log
    }

}

task bowite2_align {

    input {
        # reads pair
        File fastq1
        File fastq2

        # pre-built genome index files
        Array[File] genome_indexes = 
            ["oss://genean-seqmeta/pipeline/hg38/hg38.1.bt2","oss://genean-seqmeta/pipeline/hg38/hg38.2.bt2",
            "oss://genean-seqmeta/pipeline/hg38/hg38.3.bt2","oss://genean-seqmeta/pipeline/hg38/hg38.4.bt2",
            "oss://genean-seqmeta/pipeline/hg38/hg38.rev.1.bt2","oss://genean-seqmeta/pipeline/hg38/hg38.rev.2.bt2",
            "oss://genean-seqmeta/pipeline/hg38/index.sh"]

        # Resource
        Int cpu = 32
        String memory = "32G"
        String disks = "local-disk 40 cloud_ssd"
        
        String processdir
        String sample_id

        # docker image
        String docker_img = "registry.cn-shenzhen.aliyuncs.com/oseq-rd/seqmeta:0.1"
    }
    
    
    String directory =  processdir+"1.assay/02.rmhost/"
    String log_directory = processdir+"1.assay/logs/02.rmhost/"
    String log = log_directory+sample_id+".rmhost.log"
    
    String out1_name = directory+sample_id+".rmhost.1.fq.gz"
    String out2_name = directory+sample_id+".rmhost.2.fq.gz"
    String out1_name_downsize = directory+sample_id+".downsize.rmhost.1.fq.gz"
    String out2_name_downsize = directory+sample_id+".downsize.rmhost.2.fq.gz"
    
    String out = directory+sample_id+".rmhost.all.fq.gz"

    
    command <<<
        source /root/miniconda3/bin/activate
        conda activate seqmeta
        cpu_cores=$(nproc)
        mkdir -p ~{directory}
        mkdir -p ~{log_directory}
        bowtie2 --end-to-end --very-sensitive -p 30 -x \
        ~{sub(genome_indexes[0],".1.bt2","")} -1  ~{fastq1} -2 ~{fastq2}  \
        2> ~{log} | samtools fastq -N -c 5 -f 12 -1 ~{out1_name} -2 ~{out2_name} -
        seqtk sample -s99 ~{out1_name} 20000000 |gzip > ~{out1_name_downsize}
        seqtk sample -s99 ~{out2_name} 20000000 |gzip > ~{out2_name_downsize}
        cat ~{out1_name_downsize} ~{out2_name_downsize} > ~{out}
    >>>

    runtime {
        cpu: cpu
        memory: memory
        disks: disks
        docker: docker_img
    }

    output {
        Pair[File, File] rmhost_fastq_se = (out1_name, out2_name)
        File rmhost_fastq = out
        File log = log
    }
}



task hpv {

    input {
  
        Pair [File, File] rmhost_fastq_se

        # Resource
        Int cpu = 8
        String memory = "8G"
        String disks = "local-disk 40 cloud_ssd"
        
        String processdir
        String sample_id

        # docker image
        String docker_img = "registry.cn-shenzhen.aliyuncs.com/oseq-rd/seqmeta:0.5"
        File HPViewer_master = "oss://genean-seqmeta/pipeline/HPViewer-master/"
    }
    
    String directory =  processdir+"1.assay/04.hpv/"
    String hpv_out_directory = directory+sample_id
    String hpv_out = hpv_out_directory+"/"+sample_id+"_HPV_profile.txt"
    String log_directory = processdir+"1.assay/logs/02.rmhost/"
    String seqkit_result = log_directory+sample_id+".rmhost.reads.summary"

    

    command <<<
        source /root/miniconda3/bin/activate
        conda activate seqmeta
        mkdir -p ~{log_directory}
        mkdir -p ~{directory}
        seqkit stat ~{rmhost_fastq_se.left} ~{rmhost_fastq_se.right}  > ~{seqkit_result}
        python ~{HPViewer_master}/HPViewer.py -p 8 \
        -1 ~{rmhost_fastq_se.left} -2 ~{rmhost_fastq_se.right} -o ~{hpv_out_directory} \
        -d /root/miniconda3/envs/seqmeta/bin/
    >>>

    runtime {
        cpu: cpu
        memory: memory
        disks: disks
        docker: docker_img
    }

    output {
        File seqkit_result = seqkit_result
        File hpv_out = hpv_out
        File summary_L1 = directory+sample_id+"/temp/"+sample_id+"_summary_L1.txt"  
    }
}


task humann2 {

    input {
        File rmhost_fastq
  
        # Resource
        Int cpu = 8
        String memory = "16G"
        String disks = "local-disk 100 cloud_ssd"
        
        String processdir
        String sample_id

        # docker image
        String docker_img = "registry.cn-shenzhen.aliyuncs.com/oseq-rd/seqmeta:0.3"
        
        File protein_database = "oss://genean-seqmeta/pipeline/uniref90/"
        File nucleotide_database = "oss://genean-seqmeta/pipeline/chocophlan/"
        File bowtie2db = "oss://genean-seqmeta/pipeline/db_v20/"
        
    }
    String directory =  processdir+"1.assay/03.humann2/"+sample_id
    String metaphlan_bowtie2 = directory+"/"+sample_id+".rmhost.all_humann2_temp/"+sample_id+".rmhost.all_metaphlan_bowtie2.txt"
    String metaphlan_bugs = directory+"/"+sample_id+".rmhost.all_humann2_temp/"+sample_id+".rmhost.all_metaphlan_bugs_list.tsv"
    String marker_ab = directory+".marker_ab.tsv"

    command <<<
        source /root/miniconda3/bin/activate
        conda activate seqmeta
        cpu_cores=$(nproc)
        mkdir -p ~{directory}
        humann2 -i ~{rmhost_fastq} -o ~{directory} \
            --nucleotide-database ~{nucleotide_database} \
            --protein-database ~{protein_database} \
            --thread $cpu_cores \
            --metaphlan /root/miniconda3/envs/seqmeta/bin/ \
            --bowtie2 /root/miniconda3/envs/seqmeta/bin//bowtie2 \
            --metaphlan-options "--bowtie2db ~{bowtie2db}"
        python /root/miniconda3/envs/seqmeta/bin/metaphlan2.py \
            ~{metaphlan_bowtie2} --input_type bowtie2out \
            -o  ~{marker_ab} -t marker_ab_table  \
            --bowtie2db ~{bowtie2db}
    >>>

    runtime {
        cpu: cpu
        memory: memory
        disks: disks
        docker: docker_img
    }

    output {
        Array[File] humann2_result = [directory+"/"+sample_id+".rmhost.all_genefamilies.tsv", directory+"/"+sample_id+".rmhost.all_pathabundance.tsv", directory+"/"+sample_id+".rmhost.all_pathcoverage.tsv"]
        File metaphlan_bugs = metaphlan_bugs
        File metaphlan_bowtie2 = metaphlan_bowtie2
        File marker_ab = marker_ab
        File log = directory+"/"+sample_id+".rmhost.all_humann2_temp/"+sample_id+".rmhost.all.log"
    }
}