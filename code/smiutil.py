from rdkit.Chem import MolFromSmiles as smi2mol
from rdkit.Chem.rdMolDescriptors import  GetHashedMorganFingerprint,GetMorganFingerprintAsBitVect
from rdkit.DataStructs import ConvertToNumpyArray

import re

def canon_smiles(smi):
    try:
        m = smi2mol(smi)
    except:
        m = False
        print('hola cannon'+smi)
        
    if m is False:
        return False
    else:
        try:
            sim = Chem.MolToSmiles(m, isomericSmiles=True, canonical=True)
        except:
            sim = False
            #print('hola cannon2'+smi)
        return sim
    

def CleanSMI(smi):
        try: 
            clean=re.sub(r'[<>%\\/?\|]+', '', smi)
        except:
            #print('holaclean ' +smi )
            clean =False
        return clean

def HardValidSMI(smi):
    """
    A rule based function to validate a given smile string. 
    Return type: Boolean
    True: If a match is found. 
    False: Charges, Ions and No Conjugated regions found.
    """

    mysmile = CleanSMI(smi)

    if mysmile is not False:

        try: 
            illegalstring = re.search(r'\\|/|\*|Fe|\+\+|\.|\|',mysmile) #--> Sanity check!
        except:
            #print('holaill ' +smi )
            illegalstring = True
    else:
        return False
    
    if illegalstring:
        return False
    else:
        cansmile = canon_smiles(mysmile)
        if cansmile is False:
            return False
        match = re.search(r'\[\w{1,3}[\+-\.\d]+\]|\[\w{2}\]|\.|\(\*\)',cansmile)
        if match:
            return False
        else:
            conjuated = re.search(r'[a-z\W]\d+[a-zD-Z\W]+\d',mysmile) #r'[a-z\W]\d+[\w\W]+\d' or use (r'[a-z\W]\d+[a-zD-Z\W]+\d',mysmile) )
            if conjuated:
                return True
            else:
                return False


def applyMorganFP(m,**kwargs):
    fptype='bit'
    
    if 'fptype' in kwargs:
        fptype=kwargs['fptype']
    if 'fp_args' in kwargs:
        fp_args=kwargs['fp_args']     
    #fp_args = self.meta_data['fp_args']
    #fptype = self.meta_data['fptype']
    arr = np.zeros((1,))
    if fptype == 'bit': 
        arr = np.zeros((1,))
        #ConvertToNumpyArray(GetHashedMorganFingerprint(m, **fp_args), arr)
        try:
            arr = np.array(GetMorganFingerprintAsBitVect(m, **fp_args))
        except:
            print(Chem.MolToSmiles(m))
    elif fptype == 'count':
        #arr = np.zeros((1,))
        ConvertToNumpyArray(GetHashedMorganFingerprint(m, **fp_args), arr)
    return arr