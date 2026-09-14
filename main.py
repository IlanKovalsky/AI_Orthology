import argparse
import sys

from TreeNode import *
from database import *

id_counter: int = 0
BUFFER_SIZE = 100000


def generate_random_tree(n:int =1) -> TreeNode:
    """ Generates random binary tree structure\n
        input: n - number of leaves\n
        output: root of the tree """
    node:TreeNode = TreeNode(Genome())
    if n == 1:
        global id_counter
        node.leaf_id =id_counter
        id_counter+= 1
        return node
    leftLeaves:int = random.randint(1,n - 1)
    rightLeaves:int = n - leftLeaves
    node.left = generate_random_tree(leftLeaves)
    node.right = generate_random_tree(rightLeaves)
    return node


def generate_edge_lengths_exponential(current: TreeNode|None = None, mean: float= 0.1):
    """ Generates lengths of edges using exponential distribution
        input: current tree node and mean edge length
        output: none"""
    if current is None:
        return
    _lambda: float = 1.0 / mean
    while True:
        uniform: float = random.random()
        current.edge_length = -math.log(1 - uniform) / _lambda
        if EvolutionTree.get_probability_from_length_jukes_cantor(current.edge_length) < 0.5:
            break
    generate_edge_lengths_exponential(current.left, mean)
    generate_edge_lengths_exponential(current.right, mean)


def generate_simulation_seed() -> int:
    """Generates random seed\n
        input: none\n
        output: random seed """
    return random.randint(0,sys.maxsize)


def write_simulation_result(tree:EvolutionTree, ids: dict[str, int], seed: int):
    """ Writes results into file\n
        input: tree, ids states and seed of simulation\n
        output: none"""
    with open(SIMULATION_RESULTS_PATH, "a") as file:
        buffer = []
        for line in summarise_simulation(tree.first_cog_id, tree.ancestor_id, tree.tree_id, str(seed),
                         tree.root,0,ids,tree.root):
            buffer.append(line)
            if len(buffer) >= BUFFER_SIZE:
                file.writelines(buffer)
                buffer.clear()
        file.writelines(buffer)



def summarise_simulation(cog_id: int, ancestor_id: int, tree_id: int, seed: str,
                         current: TreeNode, distance: float,ids: dict[str, int], root:TreeNode):
    """ Concludes results of the simulation\n
        input: ids of COG, ancestor, tree, seed that was used for simulation,
        current node in the tree, distance from the root next ids and root node\n
        output: yields each line that needs to be written
        """
    if current is None:
        return
    distance += current.edge_length
    if current.left is None and current.right is None: #leaf
        genome_id = ids["next_genome_id"]
        ids["next_genome_id"] = genome_id + 1
        for i, gene in enumerate(current.val.genes):
            if i == 0:
                current.val.start_gene_id = ids["next_gene_id"]
            gene_id = ids["next_gene_id"]
            ids["next_gene_id"] = gene_id + 1
            line: str = ("gene_" + str(gene_id) + ",COG_" + str(cog_id + i) + ",G_" +
                                     str(genome_id) + "," + gene + ",ancestor_" + str(ancestor_id + i)
                                     + "," + str(distance) + ",tree_" + str(tree_id) + "," + seed +
                                     "," + root.val.genes[i] + "\n")
            yield line
    yield from summarise_simulation(cog_id, ancestor_id, tree_id, seed, current.left, distance, ids, root)
    yield from summarise_simulation(cog_id, ancestor_id, tree_id, seed, current.right, distance, ids, root)


def init_results_files():
    """ Check that all results files initialized"""
    if not is_file_exists(ID_STATES_PATH):
        update_id_state(get_default_id_state())
    if not is_file_exists(SIMULATION_RESULTS_PATH):
        init_simulation_results()
    if not is_file_exists(POSITIVE_PAIRS_RESULTS_PATH):
        init_pairs_results(POSITIVE_PAIRS_RESULTS_PATH)
    if not is_file_exists(NEGATIVE_PAIRS_RESULTS_PATH):
        init_pairs_results(NEGATIVE_PAIRS_RESULTS_PATH)


def write_pairs(tree: EvolutionTree, output_path: str, positive: bool):
    """Writes all positive pairs from current simulation\n
        input: evolution tree, output file path and are the pairs positive\n
        output: none"""
    leaves: list[TreeNode] = tree.root.get_leaves()
    if positive:
        label:str = '1'
    else:
        label:str = '0'
    with open(output_path,"at", buffering= 1024 * 1024) as file:
        buffer = []
        for i, leaf1 in enumerate(leaves):
            for j, leaf2 in enumerate(leaves[i + 1:], start=(i + 1)):
                dist:float =  tree.dist[i][j]
                genome1: int = tree.first_genome_id + i
                genome2: int = tree.first_genome_id + j
                m = len(leaf1.val.genes)

                if positive:
                    gene_pairs = ((k, k) for k in range(m))
                else:
                    gene_pairs = ((k, l) for k in range(m) for l in range(m) if k != l)

                for k,l in gene_pairs:
                    cog1: int = tree.first_cog_id + k
                    cog2: int = tree.first_cog_id + l
                    line:str = ("gene_" + str(leaf1.val.start_gene_id + k) + ",gene_"
                                + str(leaf2.val.start_gene_id + k) + ",COG_" +
                                str(cog1) + ",COG_" + str(cog2) + "," + label + "," + str(dist)
                                + ",G_" + str(genome1) + ",G_" + str(genome2) + "\n")
                    buffer.append(line)
                    if len(buffer) >= BUFFER_SIZE:
                        file.writelines(buffer)
                        buffer.clear()
                        print(str(i) + "," + str(j) + "," + str(k) + "," + str(l))
        file.writelines(buffer)


def write_positive_pairs(tree: EvolutionTree):
    """Writes all positive pairs of the simulation"""
    write_pairs(tree,POSITIVE_PAIRS_RESULTS_PATH, True)


def write_negative_pairs(tree: EvolutionTree):
    """Writes all negative pairs of the simulation"""
    write_pairs(tree, NEGATIVE_PAIRS_RESULTS_PATH, False)

def get_arguments():
    parser = argparse.ArgumentParser(description="Generate random mutation tree")
    parser.add_argument("--n", help="number of leaves", default=100, type=int)
    parser.add_argument("--m", help="number of COGs", default=1000, type=int)
    parser.add_argument("--l", help="length of gene", default=2000, type=int)
    parser.add_argument("--mean", help="mean edge length", default=0.1, type=float)
    parser.add_argument("--neg", help="if true then stores also the"
                                      " negative pairs of a simulation", default = False, type=bool)
    args = parser.parse_args()
    return args


def main():
    args = get_arguments()
    #generating structure of a tree and its edges
    root:TreeNode = generate_random_tree(args.n)
    tree: EvolutionTree = EvolutionTree(root, args.n)
    print("finished generating structure")

    generate_edge_lengths_exponential(root.left,args.mean)
    generate_edge_lengths_exponential(root.right, args.mean)
    print("finished generating lengths")

    tree.calculate_leaf_distances()

    #setting random seed
    random_seed = generate_simulation_seed()
    random.seed(random_seed)

    init_results_files()
    #setting tree id
    ids: dict[str,int] = get_id_state()
    tree.tree_id = ids["next_tree_id"]
    ids["next_tree_id"] = tree.tree_id + 1

    #beginning of simulation
    root.val.generate_random_genome(args.m,args.l)
    tree.ancestor_id = ids["next_ancestor_id"]
    ids["next_ancestor_id"] = tree.ancestor_id + args.m
    tree.first_cog_id = ids["next_cog_id"]
    ids["next_cog_id"] = tree.first_cog_id + args.m
    tree.first_gene_id = ids["next_gene_id"]
    tree.first_genome_id = ids["next_genome_id"]

    tree.mutate_tree_jukes_cantor()
    print("finished mutating tree")

    write_simulation_result(tree, ids, random_seed)
    print("finished writing results in table")
    #print(root.convert_to_newick_format(True) + ";")
    #tree.print_tree_dfs()
    update_id_state(ids)

    write_positive_pairs(tree)
    print("finished writing positive pairs")
    if args.neg:
        write_negative_pairs(tree)
    tree.write_tree_newick_format("tree_" + str(tree.tree_id) + ".tree")



if __name__ == "__main__":
    main()