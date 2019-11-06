##
setClass("Person",representation(name="character",
                                    age="numeric"))
setClass("Employee",representation(boss="Person"),
contains="Person")

hadley <- new("Person",name="Hadley",age=31)

hadley <- new("Person",name="Hadley",age="thirty")

hadley <- new("Person",name="Hadley",sex="male")

hadley <- new("Person",name="Hadley")
hadley@age
slot(hadley,"age")

setClass("Person",representation(name="character",
                                 age="numeric"),
         prototype(name=NA_character_,age=NA_real_))
hadley <- new("Person",name="Hadley")
hadley@age

getSlots("Person")

check_person <- function(object){
  errors <- character()
  length_age <- length(object@age)
  if(length_age != 1){
    msg <- paste("Age is length",length_age,". Should be 1", sep="")
    errors <- c(errors,msg)
  }

  length_name <- length(object@name)
  if(length_name != 1){
    msg <- paste("Name is length ",length_name,". Should be 1", sep="")
    errors <- c(errors,msg)
  }
  
  if(length(errors) == 0 )TRUE else errors
}

setClass("Person",representation(name="character",age="numeric"),
         validity=check_person)
new("Person",name="Hadley")

new("Person",name="Hadley",age=1:10)

hadley <- new("Person",name="Hadley",age=31)
hadley@age <- 1:10
#check
validObject(hadley)

sides <- function(object)0
setGeneric("sides")

setGeneric("sides",function(object){
  standardGeneric("sides")
})

setClass("Shape")
setClass("Polygon",representation(sides="integer"),contains="Shape")
setClass("Triangle",contains="Polygon")
setClass("Square",contains="Polygon")
setClass("Circle",contains="Shape")

setMethod("sides",signature(object="Polygon"),function(object){
  object@sides
})

setMethod("sides",signature("Triangle"),function(object)3)
setMethod("sides",signature("Square"),function(object)4)
setMethod("sides",signature("Circle"),function(object)Inf)

setGeneric("sides",valueClass = "numeric",function(object){
  standardGeneric("sides")
})
setMethod("sides",signature("Triangle"),function(object)"three")
sides(new("Triangle"))

showMethods("sides")
showMethods(class="Polygon")

##Method dispath
setClass("A")
setClass("A1",contains="A")
setClass("A2",contains="A1")
setClass("A3",contains="A2")

setGeneric("foo",function(a,b)standardGeneric("foo"))
setMethod("foo",signature("A1","A2"),function(a,b)"1-2")
setMethod("foo",signature("A2","A1"),function(a,b)"2-1")

foo(new("A2"),new("A2"))

setMethod("foo",signature("A2","A2"),function(a,b)"2-2")
foo(new("A2"),new("A2"))

setGeneric("type",function(object)standardGeneric("type"))
setMethod("type",signature("matrix"),function(object)"matrix")
setMethod("type",signature("character"),function(object)"character")

type(letters)
type(matrix(letters,ncol=2))

foo <- structure(list(x=1),class="foo")
type(foo)

setOldClass("foo")
setMethod("type",signature("foo"),function(x)"foo")
type(foo)

setMethod("+",signature(e1="foo",e2="numeric"),function(e1,e2){
  structure(list(x=e1$x+e2),class="foo")
})
foo+3

##Inheritance
setClass("Vehicle")
setClass("Truck",contains="Vehicle")
setClass("Car",contains="Vehicle")

setClass("Inspector",representation(name="character"))
setClass("StateInspector",contains="Inspector")

setGeneric("inspect.vehicle",function(v,i){
  standardGeneric("inspect.vehicle")
})

setMethod("inspect.vehicle",signature(v="Vehicle",i="Inspector"),
                                      function(v,i){
                                        message("Looking for rust")
                                        })
setMethod("inspect.vehicle",signature(v="Car",i="Inspector"),
         function(v,i){
           callNextMethod() #perform vehicle inspection
           message("Checking seat belts")
         })

inspect.vehicle(new("Car"),new("Inspector"))

setMethod("inspect.vehicle",
          signature(v="Truck",i="Inspector"),
          function(v,i){
            callNextMethod()
            message("Checking cargo attachments")
          })
inspect.vehicle(new("Truck"),new("Inspector"))

setMethod("inspect.vehicle", 
          signature(v = "Car", i = "StateInspector"),
          function(v, i) {
            callNextMethod()
            message("Checking insurance")
          })
inspect.vehicle(new("Car"),new("StateInspector"))

##Method dispatch2
setClass("C",contains="character")
setClass("B",contains="C")
setClass("A",contains="B")

a <- new("A","a")
b <- new("B","b")
c <- new("C","c")

setGeneric("f",function(x,y)standardGeneric("f"))

setMethod("f",signature("C","C"),function(x,y)"c-c")
setMethod("f",signature("A","A"),function(x,y)"a-a")

f(c,c)
f(a,a)

f(b,b)
f(a,c)

setClass("BC",contains=c("B","C"))
bc <- new("BC","bc")

setMethod("f",signature("B","C"),function(x,y)"b-c")
setMethod("f",signature("C","B"),function(x,y)"c-b")
f(b,b)

setMethod("f",signature("C","ANY"),function(x,y)"C-*")
setMethod("f",signature("C","missing"),function(x,y)"C-?")

setClass("D",contains="character")
d <- new("D","d")

f(c)
f(c,d)

##In the wild
library(EBImage)
setClass("Image",
         representation(colormode="integer"),
         prototype(colormdoe=Grayscale),
         contains="array"
         )

imageData = function(y){
  if(is(y,'Image'))y@.Data
  else y
}

setMethod("Ops",signature(e1="Image",e2="Image"),
          function(e1,e2){
            e1@.Data=callGeneric(imageData(e1),imageData(e2))
            validObject(e1)
            return(e1)
          })
setMethod("Ops",signature(e1="Image",e2="numeric"),
          function(e1,e2){
            e1@.Data=callGeneric(imageData(e1),e2)
            validObject(e1)
            return(e1)
          })
setMethod("Ops",signature(e1="numeric",e2="Image"),
          function(e1,e2){
            e2@.Data=callGeneric(e1,imageData(e2))
            validObject(e2)
            return(e2)
          })












