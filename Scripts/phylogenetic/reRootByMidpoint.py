from ete3 import Tree
import sys
def main():

    t = Tree(str(sys.argv[1]), format=1)
    midpoint = t.get_midpoint_outgroup()
    t.set_outgroup(midpoint)
    outfile=sys.argv[1]+"_midpoint"
    t.write(format=1, outfile=str(outfile))

main()
