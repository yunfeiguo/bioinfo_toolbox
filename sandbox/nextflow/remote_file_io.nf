#!/usr/bin/env nextflow

pdb = file('http://files.rcsb.org/header/5FID.pdb')
println pdb.text
target = file('tmp.txt')
pdb.withReader { source ->
	String line
	while(line = source.readLine()) {
		target << line
	}
}
println "---------"
println target.text
