#!/usr/bin/env nextflow
square = { it * it }

println square(9)


x = [1,2,3,4].collect(square).join()

println x


printMapClosure = { key, value ->
	println "$key = $value"
}

['Yue' : 'Wu', 'Mark' : 'Williams', 'Sudha' : 'Kumari'].each(printMapClosure)

myMap = ["China": 1, "India" : 2, "USA" : 3]
result = 0
myMap.keySet().each( { result += myMap[it] } )
println result
