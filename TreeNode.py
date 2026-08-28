import random
import math

#mutation alphabet for each letter
replacement: dict[str, str] = {
    "A": "CGT",
    "C": "AGT",
    "G": "ACT",
    "T": "ACG",
}


class TreeNode:
    def __init__(self, val:Genome, left: TreeNode| None=None, right: TreeNode| None=None,
                 leaf_id: int = -1, edge_length:float = 0.0):
        self.val:Genome = val
        self.left:TreeNode| None = left
        self.right:TreeNode| None = right
        self.leaf_id:int = leaf_id
        self.edge_length: float = edge_length #from parent


    def convert_to_newick_format(self, is_root: bool = False) -> str:
        """ Writes the tree as string in Newick format\n
                input: is current node a root\n
                output: string of tree in Newick format """
        if self.leaf_id != -1:
            if is_root:
                return str(self.leaf_id)
            return str(self.leaf_id) + ":" + str(self.edge_length)
        result:str =  ("(" + self.left.convert_to_newick_format() + "," +
                self.right.convert_to_newick_format() + ")")
        if not is_root:
            result += ":" + str(self.edge_length)
        return result

    def get_leaves(self) -> list[TreeNode]:
        """ Returns list of the leaves of the tree """
        leaves: list[TreeNode] = []
        if self.left is None and self.right is None:
            leaves.append(self)
            return leaves
        leaves.extend(self.left.get_leaves())
        leaves.extend(self.right.get_leaves())
        return leaves

class Genome:
    def __init__(self, genes:list[str] |None = None):
        if genes is None:
            self.genes: list[str] = []
        else:
            self.genes = genes


    @staticmethod
    def generate_random_gene(l: int = 0) -> str:
        """ Generates a random gene with length l\n
                input: l - length of the gene\n
                output: gene as a string """
        alphabet: str = "ACGT"
        gene:str = ""
        for i in range(l):
            letter: int = random.randint(0, len(alphabet) - 1)
            gene += alphabet[letter]
        return gene


    def generate_random_genome(self, m:int = 0, l:int = 0):
        """ Generates m random genes with length l\n
                input: m - number of genes, l - length of each gene\n
                output: none """
        for i in range(m):
            self.genes.append(Genome.generate_random_gene(l))


    @staticmethod
    def mutate_gene_jukes_cantor(gene:str = "", p: float = 0) -> str:
        """ Mutates gene based on Jukes Cantor mutation process\n
                input: gene - string of the gene, p - probability of mutation of  each letter\n
                output:  string of gene after mutation """
        mutation: str = ""
        for letter in gene:
            if random.random() < p:
                alphabet: str = replacement[letter]
                index:int = random.randint(0,len(alphabet) - 1)
                mutation += alphabet[index]
            else:
                mutation += letter
        return mutation


    def mutate_jukes_cantor(self, p: float = 0):
        """ Mutates all genes of genome based on Jukes Cantor mutation process\n
                input: p - probability of mutation of  each letter\n
                output: none """
        for i in range(len(self.genes)):
            self.genes[i] = Genome.mutate_gene_jukes_cantor(self.genes[i], p)


    def __str__(self): #for debugging
        string: str = ""
        for gene in self.genes:
            string += gene + "\n"
        return string[:len(string) - 1]


class EvolutionTree:
    def __init__(self, root: TreeNode,n: int = 0, tree_id: int = -1, ancestor_id: int = -1,
                 cog_id = -1,):
        self.root: TreeNode = root
        self.tree_id: int = tree_id
        self.ancestor_id: int = ancestor_id
        self.cog_id = cog_id
        self.dist = [[0.0] * n for i in range(n)] #distance matrix


    @staticmethod
    def get_probability_from_length_jukes_cantor(length: float = 0.0):
        """ calculates probability in Jukes Cantor using length of edge
                input: length
                output: probability of mutation"""
        return 1 - (1 / 4 + 3 / 4 * math.exp(-4 * length / 3))


    def mutate_tree_jukes_cantor(self):
        """ Mutates the tree based on Jukes Cantor mutation process"""
        EvolutionTree.mutate_tree_jukes_cantor_static(self.root.left, self.root)
        EvolutionTree.mutate_tree_jukes_cantor_static(self.root.right, self.root)


    @staticmethod
    def mutate_tree_jukes_cantor_static(current: TreeNode | None = None,
                                        parent: TreeNode | None = None):
        """ Going over a tree and mutate genome at each node based on their parent genome and Jukes Cantor mutation process
                input: current node, parent of current node
                output: none"""
        if current is None:
            return
        current.val.genes = parent.val.genes.copy()
        current.val.mutate_jukes_cantor(EvolutionTree.get_probability_from_length_jukes_cantor(current.edge_length))
        EvolutionTree.mutate_tree_jukes_cantor_static(current.left, current)
        EvolutionTree.mutate_tree_jukes_cantor_static(current.right, current)


    def print_tree_dfs(self):
        """ Prints the tree in a DFS order """
        EvolutionTree.print_tree_dfs_static(self.root)


    @staticmethod
    def print_tree_dfs_static(current: TreeNode, prefix: str = "", is_last=True):
        """ Print tree based on DFS. (For tests only)
                input:  current node, prefix to current root, is it last child of parent
                output: none """
        if current is None:
            return
        print(prefix, end="")
        new_prefix: str = prefix
        if is_last:
            print("└── ", end="")
            new_prefix += "    "
        else:
            print("├── ", end="")
            new_prefix += "│   "
        print(str(current.edge_length) + "  ", end="")
        print(current.val)
        children: list = []
        if current.left is not None:
            children.append(current.left)
        if current.right is not None:
            children.append(current.right)
        for i in range(len(children)):
           EvolutionTree.print_tree_dfs_static(children[i], new_prefix, i == len(children) - 1)

    def write_tree_newick_format(self, path: str):
        """ Writes tree structure into file in Newick format
                input: root of the tree, path - where to save the file
                output:  none """
        with open(path, "wt") as file:
            file.write(self.root.convert_to_newick_format(True) + ";")


    def calculate_leaf_distances(self):
        self.calculate_leaf_distances_dfs(self.root)


    def calculate_leaf_distances_dfs(self, current: TreeNode):
        """Calculates all distances between leaves"""
        if current.left is None and current.right is None: #leaf
            return [(current.leaf_id,0.0)]

        left_leaves = []
        right_leaves = []

        if current.left is not None:
            left_leaves = self.calculate_leaf_distances_dfs(current.left)
            left_leaves = [(i, d + current.left.edge_length) for i, d in left_leaves]

        if current.right is not None:
            right_leaves = self.calculate_leaf_distances_dfs(current.right)
            right_leaves = [(i, d + current.right.edge_length) for i, d in right_leaves]

        for i, di in left_leaves:
            for j, dj in right_leaves:
                self.dist[i][j] = di + dj
                self.dist[j][i] = di + dj
        return left_leaves + right_leaves


