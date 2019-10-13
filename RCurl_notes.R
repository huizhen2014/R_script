##RCurl三大函数
##getURL()/getForm()/postForm
library(RCurl)
##getURL()
##判断url是否存在
url.exists(url="www.baidu.com")
##收集调适信息,含有三个函数，update/value/reset
##Cumulate text across callbacks (from an HTTP response)
d <- debugGatherer()
temp <- getURL(url="www.baidu.com",debugfunction=d$update,verbose=TRUE) ##设置为TRUE

names(d$value())

cat(d$value()[1]) ##服务器地址及端口号
cat(d$value()[2]) ##服务器返回头信息
cat(d$value()[3]) ##提交服务器的头信息
d$reset() ##清除d$value()
d$value() ##空信息

##查看服务器返回的头信息
##列表形式
h <- basicHeaderGatherer()
txt <- getURL(url="http://baidu.com",headerfunction=h$update)
names(h$value())
h$value()

##查看服务器返回的头信息
##字符串形式
h <- basicTextGatherer()
txt <- getURL("http://www.baidu.com",headerfunction=h$update)
names(h$value()) ##返回NULL，说明是字符串形式，没有列
h$value() ##所有内容只是一个字符串
cat(h$value()) ##cat显示，更明晰

##查看url请求的访问信息
curl <- getCurlHandle()
txt <- getURL(url="http://www.baidu.com",curl=curl)
names(getCurlInfo(curl))

getCurlInfo(curl)$response.code
getCurlInfo(curl=curl)

##设置自己的header，把系统设置为Mac OS
myheader <- c("User-Agent"="Mozilla/5.0 (iPhone; U; CPU iPhone OS 4_0_1 like Mac OS X; ja-jp) AppleWebKit/532.9 (KHTML, like Gecko) Version/4.0.5 Mobile/8A306 Safari/6531.22.7",
              "Accept"="text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
              "Accept-Language"="en-us",
              "Connection"="keep-alive",
              "Accept-Charset"="GB2312,utf-8;q=0.7,*;q=0.7"
)
d <- debugGatherer()
curl <- getCurlHandle()
tmp <- getURL(url="http://www.baidu.com",httpheader=myheader,
              debugfunction=d$update,curl=curl, verbose=TRUE)
getCurlInfo(curl)$response.code
cat(d$value()[3]) ##提交服务器的头信息

##getForm()
##bingyin内搜素'RCUL'的url为
url <- c("http://cn.bing.com/search?q=Rcurl&go=Submit&qs=n&
         form=QBLH&sp=-1&pq=rcurl&sc=8-5&sk=&
         cvid=2A9A7E1B9D474402B46819841B42FD11")
##这里q=Rcul 这里就是关键字rcurl,以？开始，后续以&为分隔符
getFormParams(query=url) ##查看url的结构和值
names(getFormParams(query=url))

##tmp <- strsplit(url,split="?",fixed=T)[[1]][2]
##value <- c()
##for(i in strsplit(tmp,"&",fixed=T)[[1]]){
##  j <- strsplit(i,"=",fixed=T)[[1]]
##  k <- paste0(j[1],"=","\"",j[2],"\"")
##  value <- append(value,k)
##}
##value=paste0(value,collapse = ",")
params <- paste0(names(getFormParams(query=url)),"=","\"",getFormParams(query=url),"\"",collapse = ",")
tmp <- getForm(uri="http://cn.bing.com/search",
               q="电影团购",go="Submit",qs="n",form="QBLH",
               sp="-1",pq="rcurl",sc="8-5",sk="",
               cvid="2A9A7E1B9D474402B46819841B42FD11")
write.table(tmp,"tmp.txt")  ##读到本地文件

##postForm()
##以保密的形式上传我们所要页面提交的信息，
##然后获取服务器端返回该页面信息。例如登陆一个页面，
##需要账户和密码，那么我们需要提交账户和密码，提交的信息要加密，
##然后抓取登陆后的页面信息。

##getBinaryURL()  下载一个文件
url <- "http://rfunction.com/code/1201/120103.R"
tmp <- getBinaryURL(url)
note <- file("120103.R",open="web")
writeBin(tmp,note)
close(note)

##getBinaryURL() 批量下载文件
url <- "http://rfunction.com/code/1202/"
tmp <- getURL(url,verbose=TRUE) ##获取网页,可能需要多试几次
##查看网页源码，之后确定抓取信息的"代码字串"特征
tmp_files <- strsplit(tmp,split="<td><a href=\"")[[1]]
tmp_files1 <- strsplit(tmp_files,split="\"")
tmp_files2 <- lapply(tmp_files1,function(x)x[1])
files <- unlist(tmp_files2)
files <- files[c(-1,-2)]

##setwd("/destnation/dir/for/download)
i=1
base=url
for(i in 1:length(files)){
  url=paste(base,files[i],sep="") ##拼接url地址
  temp <- getBinaryURL(url) ##获取网页内容
  note=file(paste("1202",files[i],sep="."),open="wb") ##文件属性
  writeBin(temp,note) ##文件写入内容
  close(note)  ##关闭文件
  Sys.sleep(2) ##防止频繁访问被拉黑，采用sleep()间隔访问
}

##XML
library(XML)
##中文界面抓出来是乱码
url <- "http://data.earthquake.cn/datashare/datashare_more_quickdata_new.jsp"
##英文界面没问题
url <- "http://219.143.71.11/wdc4seis@bj/earthquakes/csn_quakes_p001.jsp"

wp <- getURL(url)
doc <- htmlParse(wp,asText=TRUE) ##这里记得encoding
tables <- readHTMLTable(doc,header=F,which=2) ##选取第二个表格

head(tables)





