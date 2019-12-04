##R advanced / 2019.10.17

##part I/Foundations
##Chapter2/Data structures

##homogeneous/ atomic vector; matrix; array
##heterogeneous/ list; data frame
##vector/ typeof();lenght();attributes();names(); as.*; is.*
numeric <- c(1,2.5,4.5)
##L 后缀区分整数和浮点数
integer <- c(1L,6L,10L)
##使用T/F或TRUE/FALSE构建逻辑向量
logical <- c(T,FALSE,TRUE,FALSE)
character <- c("these are","some strings")
##合并向量
c(1, c(2,c(3,4)))

##检测类型
typeof(numeric)
typeof(integer)
is.integer(integer)
is.double(integer)
is.numeric(integer)
##逻辑值TRUE为1，FALSE为0
as.numeric(c(F,F,T))
##TRUE总数
sum(mtcars$cyl == 4)
##符合对象(1/0)的均值
mean(mtcars$cyl == 4)
##list 不同于c()，可包含不同类型
x <- list(1:3,"a",c(T,F,T),c(2.3,5.9))
str(x)
x <- list(list(list(list())))
str(x)
is.recursive(x) ##TRUE

##c将多个list合并为1个list
x <- list(list(1,2),c(3,4))
y <- c(list(1,2),c(3,4))
str(x)
str(y)
typeof(x)
typeof(y)
##list也用来家线性模型
mod <- lm(mpg ~ wt, data=mtcars)
is.list(mod)
names(mod)
##使用unnlist将list变为向量，类似c()
typeof(unlist(mod))
typeof(c(mod))

##属性 attr()/attributes()
y <- 1:10
attr(y,"my_atrribute") <- "This is a vector"
attr(y,"my_atrribute")
str(attributes(y))
##structure()函数返回添加属性的新对象
structure(1:10,my_attribute="This is a vector")
structure(1:6, dim = 2:3)

##names()/class()/dim()
x <- 1:3; names(x) <- c("a","b","c")
x <- setNames(1:3,c("a","b","c"))
##使用unmame()或names(x) <- NULL去除名称

##Factors
x <- factor(c("a","b","b","a"))
class(x)
levels(x)
##不能使用不在leves内的值, 同时不能合并factor
sex_char <- c("m","m","m")
sex_factor <- factor(sex_char,levels=c("m","f"))
table(sex_char)
table(sex_factor) ##同时可见f水平为0
##将数字factor恢复成数字,避免自动读取factor时非标准数字带来的问题
z <- factor(c(12,1,9))
as.numeric(z) #wrong
as.numeric(as.character(z))

##Matrices and arrays
##matrix()/array()/dim()
a <- matrix(1:6,nco=3)
b <- array(1:12,c(2,3,2))
c <- 1:6; dim(c) <- c(3,2)
##length()/nrow()/ncol()/dim()/names()/rownames()/colnames()/dimnames()
c(a,a)
cbind(a,a)
rbind(a,a)
abind:abind(b,b)
##使用dim()构建grid-like结构
l <- list(1:3,"a",T,1.0)
dim(l) <- c(2,2)
l[1,1] == l[[1]]; l[2,1] == l[[2]]

##Data frames
##2维结构，因此具有matrix和list的性质
df <- data.frame(x=1:3,y=c("a",'b','c'),stringsAsFactors = FALSE)
str(df)
typeof(df)
cbind(df,data.frame(z=3:1))
rbind(df,data.frame(x=10,y="z"))
##若合并不全都含有相同值都数据框，使用:plyr:rbind:fill()
##由于data是一向量列表，因此数据框可拥有list的列
df <- data.frame(x=1:3)
df$y <- list(1:2,1:3,1:4)
df
##若将list给到data.frame()，会将每一个term放入本身列中：
##data.frame(x=1:3,y=list(1:2,1:3,1:4)) erroe
##使用I()将list作为一个单元应用到data.frame()中
##I()，AsIs 
dfl <- data.frame(x=1:3,y=I(list(1:2,1:3,1:4)))
str(dfl)
dfl[2,"y"]
dfm <- data.frame(x=1:3,y=I(matrix(1:9,nrow=3)))
dfm

##Chpter3/Subsetting
##atomic vectors
x <- c(2.1,4.2,3.3,5.4)
x[c(3,1)]
x[c(2.1,2.9)] ##非整数截断为整数
x[-c(1,2)] ##负号表反选
x[c(T,T,F,F)] ##逻辑向量对应选择位置值
x[x>3] ##加入判断，使用逻辑对应值

##lists,同atomic vector，使用[[和$取代[

##matrices and arrays
##可指定单个向量，多个向量或矩阵来提取
##[将会提取出最低的维度单元
vals <- outer(1:5,1:5,FUN="paste",sep=",")
select <- matrix(ncol=2,byrow=TRUE,c(
  1,1,3,1,2,4
))
vals[select]

##Data frames
##同时拥有list和矩阵特征，因此提取单个向量时，
##其表现同list；提取两个向量时，表现同矩阵
df <- data.frame(x=1:3,y=3:1,z=letters[1:3])
df[df$x == 2,]

##S3对象由atomoic vector, arrays和lists组成，可使用
##以上函数和str()操作
##S4对象拥有两额外提取函数：@等同与$，slot()等同于[[

##Subsetting operators
##[,$,[[

##Simplifying vs. preserving subsetting
##simplifying [[ 返回最简单的数据结构；preserving [ 返回
##输入的数据结构的输出结果
x <- c(a=1,b=2)
str(x[1])
str(x[[1]])
z <- factor(c("a","b"))
str(z[1])
str(z[1,drop=TRUE])
a <- matrix(1:4,nrow=2)
str(a[1,,drop=FALSE])
str(a[1,])
df <- data.frame(a=1:2,b=1:2)
str(df[1])
str(df[,"a",drop=FALSE])

##$,最常见错误就是在数据框中使用变量替换：var <- "cyl",mtcars$var
##应使用mtcar[var],mtcar[[var]] == mtcars[,var,drop=T]
##但是$可以实现变量的部分匹配，可在设置中取消
##options(warnPartialMatchDollar=TRUE)
mtcars$mp
x <- list(abc=1)
x$a

##Miss/out of bounds indices
##[和[[在处理超出indices范围值是稍有差异
x <- 1:3
str(x[4])
str(x[NA_real_])
str(x[NULL])

##所有的subsetting操作都可以和分配值一起进行
##不能将整数或逻辑值和NA合并一起

##当分配值给以空的索引时，将会保存为原对象的类和结构
mtcars[] <- lappy(mtcars,as.integer) ##依然为数据框
mtcars <- lapply(mtcars,as.integer) ##变为列表

##对应list，使用NULL赋值可清楚该index;
##使用list(NULL)可仅清楚index的值，保留index

##Applications
##查找表格,字符查找
grades <- sample(3,5,rep=TRUE)
info <- data.frame(
  grade=1:3,
  desc=c("Poor","Good","Excellent"),
  fail=c(TRUE,FALSE,FALSE)
)
id <- match(grades,info$grade)
info[id,]

##Boolean algebra vs sets(logical & integer subsetting)
##which()返回TRUE的位置indeices
x1 <- 1:10 %% 2 == 0
x1 <- which(x1)
x2 <- 1:10 %% 5 == 0
x2 <- which(x2)
intersect(x1,x2)
union(x1,x2)
setdiff(x1,x2)
#xor(x1,x2) ## == setdiff(union(x1,x2),intersect(x1,x2))
setdiff(union(x1,x2),intersect(x1,x2))

##Chapter4/Vocabulary
##Chapter5/Functions
f <- function(x)x
formals(f)
body(f)
environment(f)

## Lexical scopin
##scoping is the set of rules that R applies
##to go from the symbol x, to its value 10.
##仅需要知道变量值从哪里来，无需考虑怎么来的
##先查寻当前函数，然后在查询函数定义处，依次向
##上查询当前变量值
j <- function(x){
  y <- 2
  function(){ 
    c(x, y)
  } 
}
k <- j(1) ##j(1) 无显示值
k()
##k keeps around the environment in which it was defined,
##which includes the value of y

##Functions vs. variables
##one small tweak to the rule for functions
n <- function(x)x/2
o <- function(){
  n <- 10
  n(n)
}
o()
rm(n,o)

##a fresh start
##每次函数调用的对象值都是独立的

##dynamic lookup
##通过在函数内构建全局变量可解决函数变量
##和全局变量重名问题
##1. 通过函数codetools::findGlobals()查看
##一个函数所依赖的外部函数
f <- function()x+1
codetools::findGlobals(f)
##2. 或者手动改变该函数的环境,emptyenv()
environment(f) <- emptyenv()
f()

##every operation is a function call
##everything that exists is an object
##everything that happens is a function call

##`,backtick,lets you refer to functions or variables
##that have otherwise reserved or illegal names:
x <- 10; y <- 5
x+y
`+`(x,y)
x[3]
`[`(x,3)
{print(1);print(2);print(3)}
`{`(print (1),print(2),print(3))

add <- function(x,y)x+y
sapply(1:10,add,3)
sapply(1:10,`+`,3) ##返回对象+的值
sapply(1:10,"+",3) ##根据名称"+搜索函数

##Function arguments
##calling functions
##可以通过位置，完整名称，部分名称来指定参数
##参数匹配首先根据精确名称，后根据前缀和后缀匹配。
f <- function(abcdef,bcde1,bcde2){
  list(a=abcdef,b1=bcde1,b2=bcde2)
}
str(f(3,2,a=1))
str(f(1,3,b=1))

##calling a function given a list of arguments
args <- list(1:10,na.rm=TRUE)
do.call(mean,list(1:10,na.rm=TRUE))
mean(args)

##lazy evaluation
f <- function(x){
  force(x)
  10
}
system.time(f(Sys.sleep(10)))

add <- function(x){
  #force(x)
  function(y)x+y
}
adders <- lapply(1:10,add)

f <- function(x=ls()){
  a <- 1
  x
}
f() ##evaluated inside f
f(ls()) ##evaluated in global environment

a <- NULL
if(is.null(a))stop("a is null")
!is.null(a) || stop("a is null")
!is.null(a) && a > 0

##...
f <- function(...){
  names(list(...))
}
f(a=1,b=2)
sum(1,2,na.mr=TRUE)

##special calls
##infix functions
##所有用户定义的插入式的函数必须以%开头
##R预定义的插入函数有:%%,%*%,%/%,%in%,%o%,%x%
"%+%" <- function(a,b)paste(a,b,sep="")
"new" %+% " string"
`%+%`("new"," string")

##replacement functions
##替换函数类似在其位置修饰参数,拥有特殊名称xxx<-
##常拥有两个参数(x,value)，必须返回修饰后对象
"second<-" <- function(x,value){
  x[2] <- value
  x
}
x <- 1:10
second(x)<-5L
x

"modify<-" <- function(x,position,value){
  x[position] <- value
  x
}
modify(x,1) <- 10
x ## x <- `modify<-`(x,1,10)

##return values
f1 <- function()1
f2 <- function()invisible(1)

##将一个值分配给多个变量
a <- b <- c <- d <- 2

##on.exit(),当functions终止时执行
in_dir <- function(path,code){
  cur_dir <- getwd()
  on.exit(setwd(cur_dir))
  forde(code)
}

##Chapter6
##OO field guide
##类定义了对象的行为，描述其特征，以及和其他
##类的联系，对应的方法
##S3, generic-function OO
##S4，类似S3，但是更正式，拥有正式的类定义，
##描述了每个类的代表和继承，拥有特殊的帮助函数
##Reference classes，RC， 和S3/S4不同。
##base types, the internal C-level types underlie the other OO systems

##function's type is "closure", but the primitive function's type is builtin
##f <- function(){} ; typeof(f); typeof(sum)
##查看类
library(pryr)
df <- data.frame(x=1:10,y=letters[1:10])
otype(df)
otype(df$x)
otype(df$y)

methods("mean") ##查看S3对象的方法函数
methods(class="ts") ##列出所有拥有该method的指定类
ftype(mean) ##from pryr 

##S3, methods are associated with functions, generic functions/generics
##defining classes and creating objects
foo <- structure(list(),class="foo") ##during the creation with structure()
foo <- list()
class(foo) <- "foo" ##modifying an existing objects
class(foo)
inherits(foo,"foo")

##S3 classes provide a constructor function
foo <- function(x){
  if(!is.numeric(x))stop("X must be numeric")
  structure(list(x),class="foo")
}

##S3 has no checks for correctness, you can change the 
##existing objects
mod <- lm(log(mpg) ~ log(disp),data=mtcars)
class(mod)
class(mod) <- "table"
print(mod)

##creating new methods and generics
##使用函数UseMethod(), 参数1:generic function名称
##参数2: method dispatch/or it will dispath on the first argument
##to the function
f <- function(x)UseMethod("f") ##创建新的generics
f.a <- function(x)"Class a" ##增加方法/function,类似当当标签为a时的
##函数
a <- structure(list(),class="a") ##创建对象a,定义类 class=a
class(a)
f(a) ##返回f.a， dispatch 

##向现有generic增加method
##a <- stucture(list(),class="a")
mean.a <- function(x)"a"
mean(a)

##method dispath
##UseMethod()创建函数名称向量,然后依次查询返回
##"defult" class makes it possible to set up 
##a fall back method for otherwise unknown classes.
f <- function(x)UseMethod("f")
f.a <- function(x)"Class a"
f.default <- function(x)"Unknown class"
f(structure(list(),class="a"))
f(structure(list(),class=c("b","a")))
f(structure(list(),class="c"))

##S4 works in a similar way to S3, but it adds formality 
##and rigour. Methods still belong to functions, not classes,
##but: 
##1. classes hava a formal definition, describing their fields
##and inheritance structure(parent classes)
##2. Method dispath can be based on multiple arguments to 
##generic function not just one
##There is a special operator, @ , for extracting fileds(aka slots)
##out of an S4 object.

##There aren't any S4 classes in the commonly used base packages(
##stats,graphics,utils,datasets,and base)
##Create an S4 object from the built-in stats packages
library(stats4)
library(pryr)
y <- c(26,17,13,12,20,5,9,8,5,4,8)
nLL <- function(lambda) - sum(dpois(y,lambda,log=TRUE))
fit <- mle(nLL, start=list(lambda=5),nobs=length(y))
isS4(fit)
otype(fit)
isS4(nobs)
ftype(nobs)

foo <- structure(list(),class="foo")
foo <- list()
class(foo) <- "foo"
inherits(foo,"foo")

f <- function(x)UseMethod("f")
f.a <- function(x)"Class a"
a <- structure(list(),class="a")
class(a)
f(a)

mean.a <- function(x)"a"
mean(a)


f <- function(x)UseMethod("f")
f.a <- function(x)"Class a"
f.default <- function(x)"Unknown class"

f(structure(list(),class="a"))
f(structure(list(),class=c("b","a")))
f(structure(list(),class="c"))
  
setClass("Person",
         slots=list(name="character",age="numeric"))
setClass("Employee",
         slots=list(boss="Person"),contains="Person")
alice <- new('Person',name="Alice",age=40)
john <- new("Employee",name="John",age=20,boss=alice)

alice@age
slot(john,"boss")

setClass("RangedNumeric",
         contains="numeric",
         slots=list(min="numeric",max="numeric"))
rn <- new("RangedNumeric",1:10,min=1,max=10)
rn@min
rn@.Data

setGeneric("union")
setMethod("union",c(x="data.frame",y="data.frame"),
          function(x,y){
            unique(rind(x,y))
          })

setGeneric("myGeneric",function(x){
  standardGeneric("myGeneric")
})

selectMethod("nobs",list("mle"))
pryr::method_from_call(nobs(fit))

##Environment
e <- new.env()
parent.env(e)
identical(e, globalenv())
ls(e)

ls(e)
e$a <- 1
ls(e)
e$a
e$.b <- 2
ls(e)
ls(e, all=TRUE)
as.list(e)
str(as.list(e,all=TRUE))

e$b <-2 
e$b
e[["b"]]
get("b",e)

x <- 1
e1 <- new.env()
get("x",e1)

e2 <- new.env(parent=emptyenv())
get("x",e2)

exists("b",e)
exists("b",e,inherits = FALSE)
exists("b",e,inherits = FALSE)

globalenv()
baseenv()
emptyenv()

search()

as.environment("package:stats")
ls()

library(pryr)
where("where")
where("mean")
where("t.test")

F <- function(env=parent.frame()){
  if(identical(env,emptyenv())){
    return(emptyenv())
    }else{
      tmp=1
      }
  if(tmp){
    print(env)
    F(env=parent.env(env))
  }
  }
F()

where("mean")
W <- function(x){
  tmp <- as.environment(where("x",env=parent.env(where("x"))))
  cat(x,"\t")
  tmp
  
}

##3
G <- function(name,env){
  if(is.na(match(name,ls(env)))){
    return("Break")
  }else if(! is.na(match(name,ls(env)))){
    if(! is.na(match(name,search()))){
      print("good")
      
    }
  }
  }
  
  if(match(name,ls()) && !(match(name,fun))){
    G(name,env=parent.env(env))
  }else if(match(name,ls()) && match(name,fun)){
    return(name,match(name,ls(),name,fun))
  }

  
if( is.na(match("base",search())) ){
  c <- ls(paste0("package",":","base"))
  }

library(pryr)
y <- 1
f <- function(x)x+y
environment(f)
where("f")
environment(plot)
environment(t.test)
where("plot")

environment("a")
environment(x)

plus <- function(x){
  function(y)x+y
}
plus_one <- plus(1)
plus_one(10)
plus_two <- plus(2)
plus_two(10)

funenv <- function(f){
  f <- match.fun(f)
  environment(f)
}

f <- function(x)x+y
environment(f)

f <- function(){
  x <- 10
  function(){
    x}
}

f2 <- function(){
  x <- 10
  function(){
    def <- get("x",environment())
    cll <- get("x",parent.frame())
    list(defined=def,called=cll)
    }
}
g2 <- f2()
#x <- 20
str(g2())

plus <- function(x){
  function(y)x+y
}

plus_one <- plus(1)
plus_one(10)
plus_two <- plus(2)
plus_two(10)
environment(plus)
environment(plus_one)
parent.env(environment(plus_one))

f <- function(x)x+y
environment(f) <- emptyenv()
f(1)

f <- function(){
  list(
    e <- environment(),
    p <- parent.env(environment())
  )
}
str(f())
str(f())
funenv("f")

f <- function(){
  x <- 10
  function(){
    x
  }
}
g <- f()
x <- 20
g()

##local
df <- local({
  x <- 1:10
  y <- runif(10)
  data.frame(x=x,y=y)
})
##equivalent to 
df <- (function(){
  x <- 1:10
  y <- runif(10)
  data.frame(x=x,y=y)
})()

a <- 10
my_get <- NULL
my_set <- NULL
local({
  a <- 1
  my_get <<- function()a
  my_set <<- function(value)a <<- value
})
my_get()
my_set(13)
my_get()

##assignment
`  ` <- "space"
`a+b` <- 3
`:)` <- "simle"
ls()

##constants
x <- 10
lockBinding(as.name("x"),globalenv())
x <- 15
rm(x)

x %<c-% 20
x <- 30
rm(x)


