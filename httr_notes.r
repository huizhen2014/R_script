##背景知识
##https://code.tutsplus.com/tutorials/http-the-protocol-every-web-developer-must-know-part-1--net-31177
##https://www.jmarshall.com/easy/http/

##httr
##针对http最重要的两部分：request，向服务器发送数据；response，服务器返回的数据

##1. httr basics
##首先家在httr，然后使用GET()函数,url
library(httr)
r <- GET("http://httpbin.org/get")
##这里就给定一个相应对象，打印该相应对应将给出一些有用信息：
##the actual url used; the http status; the file(content)type, the size
##and if it's a text file, the first few lines of output
r
##可使用多种帮助方法获取相应对象的重要部分，或直接查询对象
status_code(r)
headers(r)
str(content(r))
##同时可以使用HEAD()/POST()/PATCH()/PUT()/DELETE()等函数
##可能做习惯与GET()/POST()，前者当请求一个网页是，浏览器使用
##后者常用于向服务器提交一个格式. PUT()/PATCH()/DELETE()常被
##web APIs使用

##2. The response
##服务器返回的数据包含3部分；the status line; the header/the body
##最重要的status line部分为http status code，它告诉你请求是否成功

##The status code
##由三个数字组成，告知请求是否成功
r <- GET("http://httpbin.org/get")
http_status(r)
##成功的网页请求一般返回状态200，普遍错误为404(file not found)
##或403(permission denied)。若请求web APIs，可能返回500，为一般
##失败code(and thus not very helpful)
##the most memorable guides http status cats:
##https://www.flickr.com/photos/girliemac/sets/72157628409467125
##同时可设置当请求不成功是返回警告或错误,强烈建议使用任一，以便
##及时发现错误
warn_for_status(r)
stop_for_status(r)

##The body
##三种方法访问请求的主体内容，都采用content()函数
##content(r, "text"), 以字符向量形式访问主体
r <- GET("http://httpbin.org/get")
content(r, "text")
##httr自动使用content-type HTTP header提供的encoding解码来自服务器的content
##但是不能总是相信服务器提供信息，可指定encoding类型
content(r, "text", encoding = "ISO-8859-1")
##针对非text请求，使用raw vector访问请求主体
content(r, "raw")
##This is exactly the sequence of bytes that the web server sent, 
##so this is the highest fidelity way of saving files to disk
bin <- content(r, "raw")
writeBin(bin,"myfile.txt")

##httr针对一般的文件格式提供一系列默认的解析方式
##JSON 自动解析为命名的列表
str(content(r, "parsed"))
##The headers,http headers 对大小写不敏感
headers(r)
headers(r)$date
headers(r)$DATE

##Cookies
##使用类似方法访问cookies,储存在用户本地终端上的数据,缓存
r <- GET("http://httpbin.org/cookies/set", query=list(a=1))
cookies(r)
##Cookies are automatically persisted between requests to the same domain
r <- GET("http://httpbin.org/cookies/set", query=list(b=1))
cookies(r)

##3. The requst
##和response一样，request，请求同样包含三部分：a status line
##header/a body；status line定义了http方法(GET, POST, DELETE, etc)
##和url。同时可以向服务器发送额外数据, in the url(with query string)
##in the headers(including cookies)和in the body of POST(), PUT(),PATCH()请求

##The url query string
##向服务器发送简单对key-value对方式是query string
##例如, http://httpbin.org/get?key=val， 可使用query名称列表
##同时空值将自动从list舍弃
r <- GET("http://httpbin.org/get",
         query = list(key1="value1", "key2"="value2"))
content(r)$args

##Custom headers，使用add_headers()指定定制headers
r <- GET("http://httpbin.org/get", add_headers(Name="Hadley"))
str(content(r)$headers)

##Cookies，为简单对key-value对，如同query string。可发送指定对cookies
##使用set_cookies()
r <- GET("http://httpbin.org/cookies", set_cookies('MeWant'='cookies'))
content(r)$cookies

##Request body
##使用POST()ing，可在请求主体中包含数据.httr允许多种方式提供该数据
r <- POST("http://httpbin.org/post", body=list(a=1,b=2,c=3))
##可以使用encode函数决定该数据如何发送至服务器
url <- "http://httpbin.org/post"
body <- list(a=1,b=2,c=3)
##Form encoded
r <- POST(url, body=body, encode="form")
##Multipart encoded
r <- POST(url, body=body, encode="multipart")
##JSON encoded
r <- POST(url, body=body,encode="json")
##使用verbose查询发送至服务器对确切内容
POST(url, body=body, encode="multipart", verbose()) # the default
POST(url, body=body, encode="form", verbose())
POST(url, body=body, encode="json", b=verbose())
##PUT()/PATCH()也可包含请求主体, 且和POST()参数一致
POST(url, body=upload_file("mypath.txt"))
POST(url, body=list(x=upload_file("mypath.txt")))
##These uploads stream the data to the server: the data 
##will be loaded in R in chunks then sent to the remote server. 
##This means that you can upload files that are larger than memory.

##4. 使用httr和rvest包进行网页爬取
##4.1 使用httr downloading a website
library(httr) ##用于make HTML GET 和 POST请求
library(rvest) ##用于解析HTML 
library(tidyr) ##用于去除NA

url <- "http://www.sunat.gob.pe/cl-at-ittipcam/tcS01Alias"
website1 <- GET(url) ##同时获得httr持有对cookies
print(content(website1))

##4.2 使用rvest包对html_nodes获得网站信息
##使用html_nodes获取网页对title和tables; 同时h3用于网站title
##table用于tables
titles <- html_nodes(content(website1),"h3")
print(html_text(titles)[[1]])
tbls <- html_nodes(content(website1),"table")
print(length(tbls))
tbl2 <- html_table(tbls[[2]],fill=T)
print(tbl2)

##4.3 使用html forms下载之前月份数据
##fThe form method used is POST . 
##We will prepare a query with the fields of the form 
##and submit that info with POST function from httr
##anho对应西班牙年/mes对应西班牙月
query <- list('mes'='10','anho'='2017')
website2 <- POST(url, body=query, encode="form")
print(content(website2))
titles <- html_nodes(content(website2),"h3")
print(html_text(titles)[[1]])
tbls <- html_nodes(content(website2),"table")
print(length(tbls))
tbl2 <- html_table(tbls[[2]],fill=T)
print(tbl2)

##4.4 将数据整理成规则的数据框格式
num.cols <- dim(tbl2)[2]
num.rows <- dim(tbl2)[1]
print(dim(tbl2))

dia <- c() ##储存天数
compra <- c() ##存储购买价格
venta <- c() ## 存储出售价格

for(i in 2:num.rows){
  for(j in 1:(num.cols/3)){
    dia <- c(dia,as.numeric(tbl2[i,(j-1)*3+1]))
    compra <- c(compra, as.numeric(tbl2[i,(j-1)*3+2]))
    venta <- c(venta, as.numeric(tbl2[i,(j-1)*3+3]))
  }
}
pen.oct.2017 <- data.frame(dia,compra,venta)
pen.oct.2017 <- pen.oct.2017 %>% drop_na() ##舍弃NA值
print(pen.oct.2017,row.names=F)
