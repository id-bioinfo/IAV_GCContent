from ete3 import Tree
import sys

def main():
    print("Tree is: "+str(sys.argv[1]))
    print("Tolerance is: "+str(sys.argv[2]))
    limit=float(sys.argv[2]) #0.00000001
    #tree = Tree("(A:0.1,B:0.000001,C:0.2,(D:0.05,E:0.0001):0.00000001);")
    treefile=sys.argv[1]
    tree = Tree(treefile, format=1)
    for node in tree.get_descendants():
        if not node.is_leaf() and node._dist <= limit:
            node.delete()
    tree.write(outfile=treefile+"_collapsed", format=1)

main()

