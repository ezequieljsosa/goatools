'''
Created on 10/12/2013

@author: Ezequiel Sosa
'''

import networkx as nx

class AnnotationTreeBuilder(object):
    '''
    Builds an AmiGO like tree, using the go terms and the gene product annotations
    * Annotations are leafs
    * Only terms with annotations in them or one of their children appears
    * The GO DAG is walked using a recursive deep first algotithm
    * When a children has more than one parent, the second time a children is procesed,
      the branch started by him is recursively copied and pasted as a new branch
     
    See http://www.sequenceontology.org/gff3.shtml especial attention on column 9
    '''
        
    COPY_CHAR = "'"
    '''
    Each time a node is replicated, this subfix is added to the id of the node. 
    Ex: A, then is copied and the copied is renamed to A'
    Then if A' is replicated, its copy is renamed to A'' and so on
    '''
        
    GFF_GO_ATRIBUTE = "GO"
    '''
    Attribute of the 9th column of a GFF3 file. The values are the GO terms that annotates that gene product
    '''    
    
    GFF_ID_ATRIBUTE = "ID" # OR Name? Name can be repeated, ID not 
    '''
    Attribute of the 9th column of a GFF3 file. ID is unique and Name is not when a feature is repeated in the genome
    ''' 
    
    
    COUNT_NODE_ATTR = "count"
    '''
    Attribute of the resulting json. Association count for a GOterm
    '''
    
    ROOT_ID = "root"
    '''
    Attribute of the resulting json. ID in the node root of the networkx graph 
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.nprocessed = {}
        '''
        It has the number of times that each node was procesed/visited, 
        the number increases when the node is in a copied branch.
        Used only in the mapping methods, to avoid the node id repetition.
        '''
        self.gene_product_by_go = {}
        '''
        Dictionary of  GOTerm > AnnotationID Qualifier List  / String => [{}1,{}2,...]
        gene_product_by_go[goterm] contains an array of each of the annotation attributes (9th column)
        '''
        self.processed = {}
        '''
        It has the annotations of the procesed nodes
        '''
        self.graph = nx.DiGraph()
        self.filtered_genes = None  
        
    def _mapping_with_increment(self,x):
        '''
        Appends prefix to all the graph nodes ids and increseases the amount of procesed times for each node
        '''
        org_label = x.replace(self.COPY_CHAR, "")
        n = self.nprocessed[org_label] if self.nprocessed.has_key(org_label)  else 0
        self.nprocessed[org_label] = n + 1         
        return self._mapping(x)
    
    def _mapping(self,x):         
        org_label = self._original_name(x)
        return org_label + (self.COPY_CHAR * self.nprocessed[org_label])
    
    def _original_name(self,label):
        '''
        Original name from a copied node
        '''
        return label.replace(self.COPY_CHAR, "")
    
    def load_annotations_by_go_dir(self,annotation_dir,go_annotated_type):
        import os
        for annotation_file in os.listdir(annotation_dir):
            self.load_annotations_by_go_file(annotation_dir + "/" +annotation_file, go_annotated_type)
    
    def load_annotations_by_go_file(self,annotation_file,go_annotated_type):
        '''
        Loads the gene_product_by_go variable by reading the GFF file and iterating over the GO attribute
        '''      
        from BCBio import GFF
        #import pprint
        #from BCBio.GFF import GFFExaminer
        # TODO poner log con datos del examiner
        #examiner = GFFExaminer()
        #in_handle = open(annotation_file)
        #pprint.pprint(examiner.parent_child_map(in_handle))
        #pprint.pprint(examiner.available_limits(in_handle))
        #in_handle.close()
        in_handle = open(annotation_file)
        limit_info = dict(gff_type = [go_annotated_type])
        for rec in GFF.parse(in_handle, limit_info=limit_info):           
            for f in rec.features:   
                if (self.GFF_GO_ATRIBUTE in f.qualifiers):
                    for q in f.qualifiers[self.GFF_GO_ATRIBUTE]:
                        if(q in self.gene_product_by_go) :
                            self.gene_product_by_go[q].append(f.qualifiers)
                        else:
                            self.gene_product_by_go[q] = [f.qualifiers]
        in_handle.close()
        
    def load_dag(self,obo_file):
        from goatools.obo_parser import GODag
        self.go_dag = GODag(obo_file)   
    
    def _clean_up_graph(self):
        #Remove nodes without annotations
        nodes_to_remove = []
        for n, d in self.graph.nodes_iter(data=True):        
            if d[self.COUNT_NODE_ATTR] == 0:
                nodes_to_remove.append(n)
                
        self.graph.remove_nodes_from(nodes_to_remove)
        
        # Remove disconected nodes                
        from itertools import chain
        #Get the list of connected components, get all but the first(biggest grapg) and then flatten the list
        disconected_nodes = list( chain(*nx.connected_components(self.graph.to_undirected())[1:]))
        # Remove the isolated nodes      
        self.graph.remove_nodes_from(disconected_nodes) 
    
    def _copy_branch(self,term,parent_id,relationship = ""):
        '''
        Copy a branch of the graph. Since all the nodes in the graph were already procesed, and its been built
        as a tree, is safe to copy its branch and append it to the parent node. 
        '''       
        from networkx.algorithms.traversal.depth_first_search import dfs_tree
        # Gets a the subtree and then a subgraph. It is used the nx.subgraph and not only dfs_tree beacause the latter does not copy the node attributes  
        subgraph = nx.subgraph(self.graph, dfs_tree(self.graph,term.id).nodes())   
        # Changes the node id of the nodes of the subgraph, so when it is joint with the main graph again, they are recognized as diferent nodes.
        subgraph = nx.relabel_nodes(subgraph,self._mapping_with_increment, copy=True) 
        # Add all the nodres from the new subgraph
        for (node, data) in subgraph.nodes_iter(data=True): 
            self.graph.add_node(node, **data)  
        #Add all the edges from the new subgraph
        self.graph.add_edges_from( subgraph.edges())
        
        branch_root = self._mapping(term.id)
        if relationship:
            self.graph.node[branch_root]["attr"] = { "rel":relationship }
        # The branch root is connected to the parent node
        self.graph.add_edge(parent_id,branch_root)
   
    def _process_relationships(self,level):    
        '''
        Iterates over all nodes and add the relationships other than is_a
        Since all the GO terms were visited, we have to copy/paste branches only
        '''
        for node in self.graph.nodes(): #For each node of the graph
            if node in self.go_dag: # if it is a valid GOTerm
                term  = self.go_dag[ self._original_name( node )]   # to obtain the term we have to use the original id                 
                if term.relationships:       
                    for related_term in [x[1] for x in term.relationships]: 
                            # Iterate over the relationship goterms. The relationship var is a duple [relationship,goterm] / <Str,Str>    
                        relationship = x[0]            
                        if related_term.id in self.processed:
                            related_term_annotations = self.processed[related_term.id]    
                            if related_term_annotations:                
                                self._copy_branch(related_term,term.id,relationship)       
                                           
                        self.processed[term.id ] = list(set(self.processed[term.id ]  + related_term_annotations)) # Remove repeated
                    self.graph.node[term.id][self.COUNT_NODE_ATTR] = len(self.processed[term.id ]) 
    
    def _order_json(self,dictionary):
        if(dictionary.has_key("children")):
            dictionary["children"] = sorted(dictionary["children"], key=lambda(node):  ( node["data"] ), reverse=True )
            for children in dictionary["children"]:
                self._order_json(children)
            
    
    def get_json(self):
        '''
        Creates the json in a tree format
        '''
        from networkx.readwrite import json_graph
        tree = json_graph.tree_data(self.graph,root="root")
        self._order_json(tree)
        return json_graph.dumps(tree)
    
    def write_file(self,file_name = 'tree.json'):
        f = open(file_name, 'w+')
        f.write(self.get_json())
        f.close() 
    

    def _process_children(self, term, level):
        '''
        Recursively process the children and get their annotations.
        '''
        annotations = []
        
        if term.children:
            for children in set(term.children):
                if children.id in self.processed: # if the child node was visited, then the branch is copied, else is recursively processed                    
                    children_annotations = self.processed[children.id]
                    if children_annotations:
                        self._copy_branch(children, term.id)
                else:
                    children_annotations = self._process_branch(children, level + 1)
                    if children_annotations:
                        self.graph.add_edge(term.id, children.id)
                annotations = annotations + children_annotations
        
        return list(set(annotations)) # Remove repeated


    def _add_annotations_to_term_node(self, term, child_annotations):
        '''
        Adds the annotations as leaf nodes, children of the current term.
        It only adds as children the annotations directly connected, not the ones in the children terms
        '''
        term_record = self.gene_product_by_go[term.id] if self.gene_product_by_go.has_key(term.id) else []
        term_annotations = []
        if term.id in self.gene_product_by_go:
            for annotation in term_record:          
                annotation_name = annotation[self.GFF_ID_ATRIBUTE][0]
                if annotation_name in child_annotations:
                    continue;  # The annotations must not be in the children       
                
                gene_symbol = annotation["gene_symbol"][0] if annotation.get("gene_symbol") else ""
                label = "-" + gene_symbol + " : " + annotation["Note"][0] # Label to be shown in the tree
                
                if not self.filtered_genes or annotation_name in self.filtered_genes:
                    # if there is a list of filtered genes, ignore the ones that are not in the list
                    annotation_node_id = self._mapping_with_increment(annotation[self.GFF_ID_ATRIBUTE][0])
                    self.graph.add_node(annotation_node_id, data=label, attr={"rel":"annotation"}, count=1)
                    self.graph.add_edge(term.id, annotation_node_id)
                    term_annotations.append(annotation_name)                                        
        
        return list(set(term_annotations).union(set(child_annotations) ))

    def _process_branch(self,term,level):     
        '''
        Creates subgraph of the term recursively
        '''    
        if not term: return
        self.graph.add_node(term.id, count=0)         
        
        # Recursive Process the children. The annotations are saved since branches with no annotations are not saved into the graph
        children_annotations = self._process_children(term, level)          
        # Adds the annotations direcly connected as graph children   
        annotations = self._add_annotations_to_term_node(term, children_annotations)   
        
        # Adds the term only if it has annnotations
        annotation_count = len(annotations)
        if annotation_count > 0:                 
            node = self.graph.node[term.id]
            node['data']= term.id + " : " + term.name + " (" +  str( annotation_count) + ")"
            node[self.COUNT_NODE_ATTR]= annotation_count    
              
        self.processed[term.id ] = annotations # Adds the term as procesed
        return  annotations
    

    def remove_count_attr(self):
        nodes_to_remove = []
        for n, d in self.graph.nodes_iter(data=True):
            del d[self.COUNT_NODE_ATTR]
        
        self.graph.remove_nodes_from(nodes_to_remove)

    def process_annotations(self):   
        '''
        Builds the entire tree, using the go terms, the gene product annotations and the filtered genes
        '''  
        assert self.go_dag, "GO terms where not loaded"   
        assert self.gene_product_by_go, "no gene product was annotated with a GO term"
                
        #GO:0008150 : biological_process
        #GO:0005575 : cellular_component 
        #GO:0003674 : molecular_function 
        root_terms = ["GO:0008150", "GO:0005575","GO:0003674"]        
        self.graph.add_node(self.ROOT_ID,data="Root", state = "open")
        
        root_count = 0
        
        for root_term in root_terms: #Iterates over each root            
            root = self.go_dag.query_term(root_term)        
            
            self.graph.add_node( root_term, data=root_term + " " + root.name)
            self.graph.add_edge(self.ROOT_ID,root_term)
            
            annotations = self._process_branch(root,1) # Process the subtree
            
            self.graph.node[root_term][self.COUNT_NODE_ATTR] = len(annotations) 
            root_count += len(annotations)
            
        self.graph.node[self.ROOT_ID][self.COUNT_NODE_ATTR] = root_count
    
        self._clean_up_graph()
        self._process_relationships(1)       
        
        
        
       
    
    def walk_graph(self,node,current_level,end_level,histogram,add_parent_levels):
        if current_level > end_level :
            return
        
        for child in self.graph.successors(node):
                child_org = self._original_name(child) 
                if child.startswith("GO"):                    
                    if child_org not in histogram :
                        if add_parent_levels or current_level == end_level :
                            histogram[child_org] =  { "data" : self.graph.node[child_org]["data"], 
                                                     "count" : self.graph.node[child_org][self.COUNT_NODE_ATTR], "level" : current_level }                    
                    self.walk_graph(child, current_level + 1, end_level, histogram,add_parent_levels)
        
        
    def histogram(self,branch,level,add_parent_levels = True):
        histogram = {}
        if add_parent_levels:
            histogram[branch] =  {"data" : self.graph.node[branch]["data"], 
                            "count" : self.graph.node[branch][self.COUNT_NODE_ATTR], "level" : 0 }
        self.walk_graph(branch,1,level,histogram,add_parent_levels)
        return histogram;
        
        
        
    
    