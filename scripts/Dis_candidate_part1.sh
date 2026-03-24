#!/bin/bash

######################################################
# 参数传递
# 新增-p参数用于控制是否执行vcf处理步骤
process_vcf=0  # 默认不执行vcf处理
while getopts "w:d:s:t:v:b:g:o:V:W:H:T:hp" opt; do
  case $opt in
    w)
      Work_path="$OPTARG"
      ;;
    d)
      Datapath="$OPTARG"
      ;;
    s)
      Scriptpath="$OPTARG"
      ;;
    t)
      Toolpath="$OPTARG"
      ;;
    v)
      Vcf_gz="$OPTARG"
      ;;
    b)
      Basepair="$OPTARG"
      ;;
    g)
      Gene_list="$OPTARG"
      ;;
    o)
      Output_prefix="$OPTARG"
      ;;
    V)
      vcfgz="$OPTARG"
      ;;
    W)
      Width="$OPTARG"
      ;;
    H)
      Height="$OPTARG"
      ;;
    T)
      Threads="$OPTARG"
      ;;
    p)
      process_vcf=1  # 当指定-p时执行vcf处理
      ;;
    h)
      display_usage
      exit 0
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      display_usage
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      display_usage
      exit 1
      ;;
  esac
done

######################################################
# 帮助文档
display_usage() {
    echo -e "\nThis script automates F1 haplotype processing including combine F1's vcffile and haplotype draw."
    # 使用 $0 获取脚本名称，实现动态显示
    echo -e "\n\tUsage: bash $(basename "$0") -w <Work_path> -d <Datapath> -s <Scriptpath> -t <Toolpath> -v <Vcf_gz> -b <Basepair>  -g <Gene_list> -o <Output_prefix> -V <vcfgz> -W <Width> -H <Height> -T <Threads> [-p]"
    echo -e "\t-w: Workpath you wanna do"
    echo -e "\t-d: Datapath you store file"
    echo -e "\t-s: Scriptpath you use"
    echo -e "\t-t: Toopath you use"
    echo -e "\t-v: Orginal Vcf_gz you wanna handle"
    echo -e "\t-b: you wanna extract round ,for example,2000:only promoter;0:only gene;2000_0;promoter+gene"
    echo -e "\t-g: Gene_list you wanna analysis"
    echo -e "\t-o: finalvcf use to analysis outputprefix"
    echo -e "\t-V: vcf_gz use to Dis_hap analysi"
    echo -e "\t-W: haplotype draw Width"
    echo -e "\t-H: haplotype draw Height"
    echo -e "\t-T: Threads that use in Dis_F1_hap analysis "
    echo -e "\t-p: Optional, execute VCF file processing steps"
    echo -e "\t-h: Display this help message"
}


# 验证所需参数是否都已提供
if [[ ! "$Work_path" || ! "$Datapath" || ! "$Scriptpath" || ! "$Toolpath" || ! "$Vcf_gz" || ! "$Basepair" || ! "$Gene_list" || ! "$Output_prefix" || ! "$vcfgz" || ! "$Width"   || ! "$Height"  || ! "$Threads"  ]]; then
    echo -e "\tERROR: Missing required parameters." >&2
    display_usage
    exit 1
fi

#1.路径参数转换，与数据路径准备

##路径参数转换
workpath=$(realpath "$Work_path")  ##工作路径
datapath=$(realpath "$Datapath")   ##数据路径
scriptpath=$(realpath "$Scriptpath")  ##脚本路径
toolpath=$(realpath "$Toolpath")   ##工具路径
##数据路径
vcf_gz="${datapath}/${Vcf_gz}"  ##变异检测完合并后的初始vcf
allgene_list="${datapath}/${Gene_list}"
snpEff_config="${datapath}/snpEff1.config"

##创建遗传距离和单倍型分析用的目录
mkdir -p "${workpath}/DIS_${Basepair}"
##进入工作目录
cd "$workpath"

#2.处理原始vcf文件 (可选项，通过-p参数控制)
if [ $process_vcf -eq 1 ]; then
    echo "Starting VCF file processing..."

    ##################
    #2.1.去掉包含缺失基因型的位点和杂合基因型的位点
    ##第一步，去除包含缺失基因型的位点
    vcftools --gzvcf "$vcf_gz" --max-missing 1 --recode --stdout > temp_no_missing.vcf

    ##第二步，检查包含缺失基因型的位点是否全部去除
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' temp_no_missing.vcf | grep -F './.'

    ##第三步，去除包含杂合基因型的位点
    bcftools view -e 'GT="het"' temp_no_missing.vcf -o output_no_het.vcf  ##提取不含杂合位点的vcf文件
    bcftools view -i 'GT="het"' temp_no_missing.vcf -o output_het.vcf     ##提取含杂合位点的vcf文件

    ##第四步，检查包含杂合基因型的位点是否全部去除
    total_sites=$(grep -v '^#' temp_no_missing.vcf | wc -l)
    no_het_sites=$(grep -v '^#' output_no_het.vcf | wc -l)
    het_sites=$(grep -v '^#' output_het.vcf | wc -l)
    echo "Total sites in temp_no_missing.vcf: $total_sites"
    echo "Sites in output_no_het.vcf: $no_het_sites"
    echo "Sites in output_het.vcf: $het_sites"
    echo "Sum of no_het + het sites: $((no_het_sites + het_sites))"  # 应等于 total_sites

    ##第五步，去除中间文件
    rm -f temp_no_missing.vcf output_het.vcf

    #2.2.对过滤完成的vcf的文件进行变异注释
    ##第一步，构建MSU注释库
    mkdir -p "${workpath}/data/MSU"
    cp "${datapath}/sequences.fa" "${workpath}/data/MSU/sequences.fa"  ##拷贝基因组序列
    cp "${datapath}/genes.gff" "${workpath}/data/MSU/genes.gff"        ##拷贝gff文件
        cp $snpEff_config ${workpath}/
    java -jar "${toolpath}/snpEff/snpEff.jar" build -c snpEff1.config -gff3 -v MSU  ##构建MSU注释库

    ##第二步，进行变异注释
    java -jar "${toolpath}/snpEff/snpEff.jar" -c snpEff1.config -ud 5000 -csvStats test1.csv -htmlStats test1.html -o vcf MSU output_no_het.vcf > "${datapath}/${Output_prefix}.SNP_INDEL.ann.vcf"

    #2.3.对变异注释完成的vcf的文件进行bcftools索引
    bcftools view -Oz -o "${datapath}/${Output_prefix}.SNP_INDEL.ann.vcf.gz" "${datapath}/${Output_prefix}.SNP_INDEL.ann.vcf"  ##压缩成gz格式
    bcftools index -t "${datapath}/${Output_prefix}.SNP_INDEL.ann.vcf.gz"  ##给vcf构建索引

    # 清理临时文件
    rm -f output_no_het.vcf
    echo "VCF file processing completed."
    echo "${Output_prefix}.SNP_INDEL.ann.vcf.gz已生成，且索引已构建"
fi

#3.进行遗传距离计算和F1单倍型分析
Workpath=`realpath ${workpath}/DIS_${Basepair}`
cd $Workpath
# 生成作业脚本
for i in $(cat "$allgene_list"); do
    echo "bash $scriptpath/Dis_F1_hap.sh -w ${Workpath} -d ${datapath} -s ${scriptpath} -t ${toolpath} -v ${vcfgz} -b ${Basepair} -g ${i} -W ${Width} -H ${Height} "
done > "Dis_${Basepair}_job"

# 执行并行作业
ParaFly -c "Dis_${Basepair}_job" -CPU "$Threads"

#清理临时文件
#rm Dis_${Basepair}_job*
