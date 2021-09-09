from pepmass.glycomass import GlycanNode
from print_tree import print_tree

class print_glycan_struct(print_tree):
    def __init__(self, glycan, **kwargs):
        if isinstance(glycan, str):
            glycan = GlycanNode.from_str(glycan)
        
        super(print_glycan_struct, self).__init__(glycan, **kwargs)
        
    def get_children(self, node):
        return node.children or []

    def get_node_str(self, node):
        return str(node.monosaccharide)


if __name__ == '__main__':
    print_glycan_struct('(N(N(H(H(H(H)))(H(H(H(H)))))))')
    
    print_glycan_struct('(N(N(H(H)(H(H(H))(H(H(H)))))))')
    
    print_glycan_struct('(N(N(H(H(H(H)))(H(H)(H(H))))))')
    