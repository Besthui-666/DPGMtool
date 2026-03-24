#!/bin/bash

######################################################
# 参数传递
Process_select=0
Select=0
while getopts "w:d:s:b:p:i:m:y:o:PSh" opt; do
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
    b)
      Basepair="$OPTARG"
      ;;
    p)
      Phenotype="$OPTARG"
      ;;
    i)
      Pheinfo="$OPTARG"
      ;;
    m)
      Mapdisinfo="$OPTARG"
      ;;
    y)
      Year="$OPTARG"
      ;;
    o)
      Prefix_output="$OPTARG"
      ;;
    P)
      Process_select=1  # 当指定-P时执行总的遗传距离合并
      ;;
    S)
      Select=1 # 当指定-S时执行遗传距离处理和t检验
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
    echo -e "\nThis script automates dis merge and analysis of pearson."
    # 使用 $0 获取脚本名称，实现动态显示
    echo -e "\n\tUsage: bash $(basename "$0") -w <Work_path> -d <Datapath> -s <Scriptpath>  -b <Basepair> -p <Phenotype> -i <Pheinfo> -m <Mapdisinfo> -y <Year> -o <Prefix_output>  [-P] [-S]"
    echo -e "\t-w: Workpath you wanna do"
    echo -e "\t-d: Datapath you store file"
    echo -e "\t-s: Scriptpath you use"
    echo -e "\t-b: you wanna extract round ,for example,2000:only promoter;0:only gene;2000_0;promoter+gene"
    echo -e "\t-p: you wanna associate Phenotype with Dis"
    echo -e "\t-i: Pheinfo used in match phenotype"
    echo -e "\t-m: Map dis to handleDis"
    echo -e "\t-y: data when you use, for instance,24/2024"
    echo -e "\t-o: result of pearson's output prefix"
    echo -e "\t-P: Optional, execute generation of Dis merge steps"
    echo -e "\t-S: Optional, execute generation of Dis processing steps"
    echo -e "\t-h: Display this help message"
}


# 验证所需参数是否都已提供
if [[ ! "$Work_path" || ! "$Datapath" || ! "$Scriptpath"  || ! "$Basepair" || ! "$Phenotype"  || ! "$Pheinfo" || ! "$Mapdisinfo" || ! "$Prefix_output" || ! "$Year" ]]; then
    echo -e "\tERROR: Missing required parameters." >&2
    display_usage
    exit 1
fi
#1.路径参数转换，与数据路径准备

##路径参数转换
workpath=$(realpath "$Work_path")  ##工作路径
datapath=$(realpath "$Datapath")   ##数据路径
scriptpath=$(realpath "$Scriptpath")  ##脚本路径

#mapdisinfo=${Mapdisinfo}
phe=${Phenotype}
pheinfo=${Pheinfo}
mapdisinfo=${Mapdisinfo}
cd ${workpath}

#2.进行候选基因筛选结果整理步骤
if [ $Process_select -eq 1 ]; then
#2.1.整理遗传距离表格
##第一步，整理成功完成遗传距离计算和F1单倍型分析的基因list
#ll -thrl ${workpath}/DIS_${Basepair} |sed 's/ /\t/g' |awk '{print$NF}' |tail -n +2 > ${workpath}/DIS_${Basepair}_genelist
find ${workpath}/DIS_${Basepair} -maxdepth 1 -type d ! -path  ${workpath}/DIS_${Basepair} |sed 's/\//\t/g' |awk '{print$NF}' > ${workpath}/DIS_${Basepair}_genelist
genelist=${workpath}/DIS_${Basepair}_genelist
##第二步，合并所有基因的遗传距离
for i in `cat ${genelist}`;do cat ${workpath}/DIS_${Basepair}/${i}/${i}_Dis ;done |grep -v 'gene_ID' |sed 's/-nan/NA/g'> DIS_${Basepair}_temp
for i in `cat ${genelist}`;do cat ${workpath}/DIS_${Basepair}/${i}/${i}_Dis ;done  |sed 's/-nan/NA/g' |head -n1 > DIS_${Basepair}_title
cat DIS_${Basepair}_title DIS_${Basepair}_temp >DIS_${Basepair}_Dis ##合并生成所有基因亲本间的遗传距离矩阵
rm DIS_${Basepair}_title DIS_${Basepair}_temp ##去除中间文件
fi

if [ $Select -eq 1 ]; then
##第三步，对遗传距离进行处理，处理成跟表型重复数对应的状态
bash $scriptpath/map_data_converter.sh  ${mapdisinfo} DIS_${Basepair}_Dis ${Year}_DIS_${Basepair}_handle_Dis  ##假设遗传距离pos_taifengA的遗传距离是0.1，但是它对应的重复有三个，则可以处理成包含三个遗传距离相同数值的列。
#2.2.组间遗传距离t检验以及均值的比值计算
##第一步，进行组间遗传距离t检验以及均值比值计算
Rscript  $scriptpath/Dis_ttest_meanratio.R ${Year}_DIS_${Basepair}_handle_Dis ${pheinfo} ${Year}_DIS_${Basepair}_disttest_mean
##第二步，处理结果文件方便最后合并
sort -k1 ${Year}_DIS_${Basepair}_disttest_mean_gene_stats.txt > ${Year}_DIS_${Basepair}_disttest_mean_gene_handle_stats.tsv
Rscript $scriptpath/peasonr.r $workpath ${Year}_DIS_${Basepair}_handle_Dis  ${phe}   ${pheinfo} ${Prefix_output}_DIS_${Basepair}
fi
#3.遗传距离与生物量表型皮尔逊相关性分析
#Rscript $scriptpath/peasonr.r $workpath ${Year}_DIS_${Basepair}_handle_Dis  ${phe}   ${pheinfo} ${Prefix_output}_DIS_${Basepair}
Rscript $scriptpath/peasonr.r $workpath DIS_${Basepair}_Dis ${phe}   ${pheinfo} ${Prefix_output}_DIS_${Basepair}
