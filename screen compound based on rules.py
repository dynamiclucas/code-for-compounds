##screen on which has Five-membered rings
from rdkit import Chem
from rdkit.Chem import rdDepictor

# 读取SMILES
smiles = [line.split(',')[0] for line in open('smiles-small.csv').readlines()[1:]] 

# 获取含5元环的分子SMILES  
five_rings = []
for smi in smiles:
    mol = Chem.MolFromSmiles(smi)
    sssr = Chem.GetSymmSSSR(mol)
    for ring in sssr:
        if len(ring) == 5:
            five_rings.append(smi)
            break

print(five_rings)


##screen on a saturated unbranched compound ##无支链化合物yst员工
from rdkit import Chem
from rdkit.Chem import rdDepictor

smiles = [line.split(',')[0] for line in open('smiles-small.csv').readlines()[1:]]

saturated = []
for smi in smiles:
    mol = Chem.MolFromSmiles(smi)
    if mol.GetNumAtoms(onlyHeavy=True) == mol.GetLongestChain():
        saturated.append(smi)

print(saturated)


##screen on a compound containing heterocyclic rings  ##筛选含有杂环的化合物yst员工
##Function
from rdkit import Chem

def filter_hetero_rings(smiles_list):
    """
    筛选出含杂环的化合物
    """
    hetero_rings = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        sssr = Chem.GetSymmSSSR(mol)
    
        has_hetero = False
        for ring in sssr:
            if len(set(ring) - set([6])) > 0:     ##rdkit中碳原子的编号为6，如果减去碳原子后，集合长度大于零（集合非空）则证明仍含有非碳原子
                has_hetero = True
                break
        
        if has_hetero:
            hetero_rings.append(smi)
            
    return hetero_rings

# 测试调用yst员工  
smiles = [line.split(',')[0] for line in open('smiles-small.csv').readlines()[1:]]
print(filter_hetero_rings(smiles))


##封装为class
class HeteroRingFilter:
    
    def __init__(self):
        pass
        
    def filter(self, smiles_list):
        hetero_rings = []
        for smi in smiles_list:
            mol = Chem.MolFromSmiles(smi)
            sssr = Chem.GetSymmSSSR(mol)
        
            has_hetero = False
            for ring in sssr:
                if len(set(ring) - set([6])) > 0:
                    has_hetero = True
                    break
            
            if has_hetero:
                hetero_rings.append(smi)
                
        return hetero_rings

# 测试调用
f = HeteroRingFilter()
smiles = [line.split(',')[0] for line in open('smiles-small.csv').readlines()[1:]]
print(f.filter(smiles))
