# AI\_Orthology



## Files Description



### python files:

'main.py' - Main file.

'TreeNode.py' - Handles most operations on the tree and its mutations.

'database.py' - Handles operations related to saving the simulation results.



### results files:

'id\_states.json' - Persistence file for saving information about IDs related to results.

'simulation\_results.csv' - Table that saves the direct results of the simulation.

'tree\_i.tree' - Files that save the structure of simulated tree in Newick format.



## How to use:

Run the main.py file in CMD\\some file editor.

How To change arguments of the simulation:

CMD - Add after main.py argument that you want to change (--n number of leaves, --m number of COGs, --l length of gene, --mean mean edge length) and number that you want to use. For example: 'main.py --n 1500' will change number of leaves to 1500.

If you want to change it directly in code then in main.py at function 'get\_arguments' change the default value of arguments as you need.

