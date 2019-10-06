##http://www.bioinfo-scrounger.com/archives/417
##Uniprot Unversal Protein 整合了swissprot/trembl/pri-psd三大数据库，
##是目前使用非常广泛的数据库
##API 通过URL以及编程命令轻松访问数据库获取所需信息
##Uniprot API 获得三部分信息
##获取单个蛋白的ID（uniprot accession id）所对应的所有信息
##ID mapping 的API
##查询访问性质的API

##访问Uniprot数据小技巧，例如获取Q9UM73蛋白的fasta序列
#http://www.uniprot.org/uniprot/Q9UM73.fasta
##若获得其全面信息，可有多种格式,txt/xml/GFF
#http://www.uniprot.org/uniprot/Q9UM73.txt
#http://www.uniprot.org/uniprot/Q9UM73.xml
#http://www.uniprot.org/uniprot/Q9UM73.GFF
##或使用URL，使用R的RCurl包，批量下载蛋白信息，避免手动去
##Uniprot网页下载

##使用Uniprot idmapping功能对应的API
##https://www.uniprot.org/help/api_idmapping
##使用R实现类似爬虫功能
##https://github.com/sehanson/Uniprot_API_R/blob/master/Uniprot_API
library(httr)
library(readr)
e_mail <- "huizhen_2014@163.com"  ##输入你的邮箱用于帮助debug
us_er <- paste0("R ",e_mail)
acc1 <- "ACC" ##输入starting accession id（例如，'ACC' for uniprot accession)
acc2 <- "KO_ID" ##输入target accession id（例如，'HPAZ-ID' for human protein atlas acccession)
fmt <- 'tab' ##输入返回格式(tab for TSV)
qry <- 'P31946,P62258' ##输入查询terms
user_agent(us_er)
r <- POST('http://www.uniprot.org/uploadlists/', 
          body = list(from= acc1, to = acc2, 
                      format = fmt, query = qry), 
          encode = "form")
#print(r)
##r变量为response格式，无法直接解读
r_content <- content(r,"text")
r_tab <- read_tsv(r_content)
##can't be parsed correctly

##programmatic access, retrieving entries via queries
##可使用query定义感兴趣的数据内容，例如所有reviewed human entries
#https://www.uniprot.org/uniprot/?query=reviewed:yes+AND+organism:9606
##或者当前uniprot/swiss-prot release：
#https://www.uniprot.org/uniprot/?query=reviewed:yes+AND+created:[current TO *]

##用于查询的URL将反馈一系列数据集名称(uniprot, uniref, uniparc, taxonomy...)
##https://www.uniprot.org/help/api_queries
##query string, query syntax, https://www.uniprot.org/help/text-search
##query fields for UniportKB...https://www.uniprot.org/help/query-fields
##format, html, tab, xls, fasta, gff, txt, xml, rdf, list, rss
##columns, 逗号分隔的列名称，查询结果显示为tabl或xls格式，https://www.uniprot.org/help/uniprotkb_column_names
##include yes/no，format设置为fasta时包含isoform 序列；format设置
##为rdf时，包含参考数据的描述
##compress yes/no，返回gzipped结果
##limit，整数，查询最大结果数目限制
##offset，整数，和limit一起使用，第一个结果的offset
##例如，查询所有人类匹配antigen词条，结果以RDF/XML和tab格式输出
##https://www.uniprot.org/uniprot/?query=organism:9606+AND+antigen&format=rdf&compress=yes
##https://www.uniprot.org/uniprot/?query=organism:9606+AND+antigen&format=tab&compress=yes&columns=id,reviewed,protein names
##例如，查询所有人类交叉PDB的词条，结果以tab格式输出，且仅包含UniprotKB和PDB识别符号
##https://www.uniprot.org/uniprot/?query=organism:9606+AND+database:pdb&format=tab&compress=yes&columns=id,database(PDB)

##列出所有人swissprot蛋白对应的KEGG id和GOid
#http://www.uniprot.org/uniprot/?query=organism:9606&format=tab&columns=id,reviewed,database(KO),go-id
#使用RCurl包抓取
library(RCurl)
library(readr)
url.exists("http://www.uniprot.org/uniprot/?query=organism:9606&format=tab&columns=id,reviewed,database(KO),go-id")
d <- debugGatherer()
tmp <- getURL("http://www.uniprot.org/uniprot/?query=organism:9606&format=tab&columns=id,reviewed,database(KO),go-id",
              debug=dverbose=TRUE)
tmp2 <- read_tsv(tmp)





















