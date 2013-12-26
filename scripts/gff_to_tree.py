#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
python %prog study.file population.file gene-association.file

This program creates a tree json representation ( AmiGO like ) from a GFF3 file with annotated GO terms.
There is a basic template with jstree in /web to test it, it takes minutes to load (the browser
 might even ask if stop the process or continue waiting, even with json trees of 50mb eventualy it was loaded )
but it loads the entire tree and sorts all the nodes

The main idea is to transform from a DAC (directed acyclic graph) to a tree, by coping branches.

The generated json is something like:

{ id: GO:0005575         <-- Node id
  data: GO:0005575 : cellular_component (37) <-- Label to show in the tree
  attr: {"rel": "is_a"}    <-- Node type / relationship with the parent, is used to select the icon
  children: [
      {
          id: Halorubrum_Sp_AJ67.transcript.471         
          data: -tbp : TATA-box-binding protein 
          attr: {"rel": "annotation"}   
      },...
  ]
} 

Parameters:
-- go_obo_file : GO term file 
-- annotation_dir: directory with all the gff3 files (with GO annotations). Must be the hierarchical files
optionals
-- annotated_type: level witch is annotated with GO terms, default: transcript
-- filtered_genes: use only the annotations that annotates these genes (uses Name gff3 attribute). Can be useful, for example, if we only want to see the tree of expressed genes

Output: tree.json: contains the json tree structure, with the GO terms as nodes and the annotations as leafs

/web/template.html: template of the tree that loads the data (ajax) of the generated tree.json. Both must be in the same web directory
"""

import os
import os.path as op
import sys
sys.path.insert(0, op.join(op.dirname(__file__), ".."))


def get_filtered_genes(genes_param):
    if "$" in genes_param:
        vec = genes_param.split("$")
        return parse_csv(vec[0],int(vec[1]),vec[2].split("="))
    else:
        return genes_param.split(",")
 
#  {"gene." => "transcript."} 
def parse_csv(csv_file,column,replacement=[]):
    import csv
    ifile  = open(csv_file, "rb")
    reader = csv.reader(ifile)
    
    filtered_genes = []

    for row in reader:
        if replacement:
            #si esta vacio
            if len(row) <= column: break;
            if not row[column].strip(): break
            vec = row[column].split(replacement[0])
            # Removes left ceros
            filtered_genes.append(vec[0] + replacement[1] +  str(int(vec[1])) )
        else:
            filtered_genes.append(row[column])
    
    ifile.close()
    return filtered_genes

if __name__ == '__main__':

    import optparse
    p = optparse.OptionParser("%prog [options] annotation_dir go_obo_file")
    p.add_option("--atype", dest="annotated_type", 
                 help="level in the gff which is annotated with GO terms, default: mRNA"
                 , action="store", type="string", default="mRNA")    
#     p.add_option("--template", dest="generate_template", action="store", type="string",
#                  default="direct", type="boolean", default=True,
#                  help=" Generates a template in generated/template.html with jstree" )
    p.add_option("--genes", dest="filtered_genes", action="store",
                 help="Only the annotations of these gene products will be used. There are 2 ways to enter this parameter:"
                 " "
                 "  **coma separated values of genes. "
                 " "
                 "  **file.csv$n$str=repl : cvs file name and nth column, from that column the gene list will be obtained."
                 " "
                 " str=repl is optional and its used when the ids from the cvs differs from the annotated ids from the gff file."
                 " For example, ids in a halorubrum.gff are 'text.transcript.othertext' and the ids in the expressed.csv are 'text.gene.othertext'"
                 " so the genes parameter is: expressed.csv$2$gene.=transcript. "
                 " Header must be removed from the csv file"
                 , type="string", default=None)
    opts, args = p.parse_args()

    # check for correct number of arguments
    if len(args) == 0:
        p.print_help()
        sys.exit(1)
    if len(args) == 1: 
        obo_file = "gene_ontology.1_2.obo"
    else:
        obo_file = args[1]
    
    assert os.path.exists(obo_file), "file %s not found!" % obo_file

    annotation_dir = args[0]
    assert os.path.exists(annotation_dir), "file %s not found!" % annotation_dir
    from goatools.annotation import AnnotationTreeBuilder
        
    builder = AnnotationTreeBuilder() 
    builder.load_dag(obo_file)
    builder.load_annotations_by_go_dir(annotation_dir,opts.annotated_type)    
    
    if opts.filtered_genes is not None:    
        builder.filtered_genes = get_filtered_genes(opts.filtered_genes)  
        
    builder.process_annotations()
    builder.remove_count_attr()
    builder.write_file()
        
    print "DONE!"
    
        
