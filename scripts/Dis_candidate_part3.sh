#!/bin/bash

######################################################
# 参数传递
while getopts "w:d:s:g:o:W:H:h" opt; do
  case $opt in
    w)
      Workpath="$OPTARG"
      ;;
    d)
      Datapath="$OPTARG"
      ;;
    s)
      Scriptpath="$OPTARG"
      ;;
    g)
      Gene_ID="$OPTARG"
      ;;
    o)
      Output_prefix="$OPTARG"
      ;;
    W)
      Width="$OPTARG"
      ;;
    H)
      Height="$OPTARG"
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
    echo -e "\nThis script automates F1 haplotype processing of posfavor."
    # 使用 $0 获取脚本名称，实现动态显示
    echo -e "\n\tUsage: bash $(basename "$0") -w <Work_path> -d <Datapath> -s <Scriptpath>  -g <Gene_ID> -o <Output_prefix> -W <Width> -H <Height> "
    echo -e "\t-w: Workpath you wanna do"
    echo -e "\t-d: Datapath you store file"
    echo -e "\t-s: Scriptpath you use"
    echo -e "\t-g: Gene_ID you wanna analysis"
    echo -e "\t-o: outputprefix"
    echo -e "\t-W: haplotype draw Width"
    echo -e "\t-H: haplotype draw Height"
    echo -e "\t-h: Display this help message"
}


# 验证所需参数是否都已提供
if [[ ! "$Workpath" || ! "$Datapath" || ! "$Scriptpath"  || ! "$Gene_ID" || ! "$Output_prefix"  || ! "$Width"   || ! "$Height" ]]; then
    echo -e "\tERROR: Missing required parameters." >&2
    display_usage
    exit 1
fi


#1.路径参数转换，与数据路径准备

##路径参数转换
workpath=$(realpath "$Workpath")  ##工作路径
datapath=$(realpath "$Datapath")   ##数据路径
scriptpath=$(realpath "$Scriptpath")  ##脚本路径

##数据路径准备
gene_ID=${Gene_ID}
phe=${datapath}/test_hap_phe.txt
subpop=${datapath}/subpop.txt

#2.提取posfavor vcf
cd ${workpath}/${gene_ID}
##第一步，生成posfavor位点文件
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' ${gene_ID}.v1.anno.vcf |sed 's/|/\//g' | awk '{d1=1;d2=1; for(i=5;i<=11;i++) if($12==$i) d1=0; for(i=13;i<=NF;i++) if($12==$i) d2=0; if(d1&&d2) print}'  |cut -f1,2 > ${gene_ID}_favorsite.txt
##判断posfavor位点文件是否为空，是空的则终止脚本
#if [ ! -s "${gene_ID}_favorsite.txt" ]; then
#        # 删除之前生成的所有相关文件
#       grep -v '#' ${gene_ID}.v1.anno.vcf |wc -l > ${gene_ID}_total_number_temp
#       cat ${gene_ID}_favorsite.txt |wc -l > ${gene_ID}_posfavor_number_temp
#       echo ${gene_ID} > ${gene_ID}_temp
#       paste ${gene_ID}_temp ${gene_ID}_total_number_temp ${gene_ID}_posfavor_number_temp >${gene_ID}_site.txt
#        rm -f ${gene_ID}_favorsite.txt ${gene_ID}*temp
#        echo "错误：上一步生成的 ${gene_ID}_favorsite.txt 为空！已删除相关文件并停止后续命令。" >&2  # 错误信息输出到 stderr
#        exit 1  # 退出脚本，状态码1（非0表示执行失败）
#fi

if [ ! -s "${gene_ID}_favorsite.txt" ]; then
    #删除之前生成的所有相关文件
    grep -v '#' ${gene_ID}.v1.anno.vcf |wc -l > ${gene_ID}_total_number_temp
    cat ${gene_ID}_favorsite.txt |wc -l > ${gene_ID}_posfavor_number_temp
    echo ${gene_ID} > ${gene_ID}_temp
    paste ${gene_ID}_temp ${gene_ID}_total_number_temp ${gene_ID}_posfavor_number_temp >${gene_ID}_site.txt
    rm -f ${gene_ID}_favorsite.txt ${gene_ID}*temp
    echo "错误：上一步生成的 ${gene_ID}_favorsite.txt 为空！已删除相关文件并停止后续命令。" >&2  # 错误信息输出到 stderr
    exit 1  # 退出脚本，状态码1（非0表示执行失败）
fi

##第二步，对${gene_ID}_newF1.vcf进行压缩和索引并基于${gene_ID}_favorsite.txt提取posfavor vcf文件
##保证pos基因型与其他都不一样--（需要考虑neu的基因型）
bcftools view -Oz -o ${gene_ID}_newF1.vcf.gz ${gene_ID}_newF1.vcf
bcftools index -t ${gene_ID}_newF1.vcf.gz
bcftools view -R ${gene_ID}_favorsite.txt ${gene_ID}_newF1.vcf.gz -o ${gene_ID}_posfavor_F1.vcf
#3.单倍型分析
#Rscript ${scriptpath}/hap.v1.r  ${workpath} ${gene_ID}_posfavor_F1.vcf ${subpop} ${phe} ${Output_prefix}_${gene_ID} ${Width} ${Height}
#3.统计总位点和posfavor位点数量
grep -v '#' ${gene_ID}_newF1.vcf |wc -l > ${gene_ID}_total_number_temp
cat ${gene_ID}_favorsite.txt |wc -l > ${gene_ID}_posfavor_number_temp
echo ${gene_ID} > ${gene_ID}_temp
paste ${gene_ID}_temp ${gene_ID}_total_number_temp ${gene_ID}_posfavor_number_temp >${gene_ID}_site.txt
rm ${gene_ID}*temp
