import sys, re, base64
import pubchempy as pcp

smiles_pattern = r'[CHONFClBrISPKacnosp0-9@+\-\[\]\(\)\/.=#$%]{5,}'

def smiles_node(name: str) -> (str):
  '''
    Queries Pubchem for the smiles string of the molecule based on the name.
      Args:
        names_list: the list of molecule names
      Returns:
        smiles_list: the list of smiles strings of the molecules    
        smiles_string: a string of the tool results
  '''
  print("smiles tool")
  print('===================================================')

  #try:
  res = pcp.get_compounds(name, "name")
  smiles = res[0].smiles
  #except:
  #  smiles = "unknown"

  return smiles


sequence_bank = {
    "MAOB": "MSNKCDVVVVGGGISGMAAAKLLHDSGLNVVVLEARDRVGGRTYTLRNQKVKYVDLGGSY\
VGPTQNRILRLAKELGLETYKVNEVERLIHHVKGKSYPFRGPFPPVWNPITYLDHNNFWR\
TMDDMGREIPSDAPWKAPLAEEWDNMTMKELLDKLCWTESAKQLATLFVNLCVTAETHEV\
SALWFLWYVKQCGGTTRIISTTNGGQERKFVGGSGQVSERIMDLLGDRVKLERPVIYIDQ\
TRENVLVETLNHEMYEAKYVISAIPPTLGMKIHFNPPLPMMRNQMITRVPLGSVIKCIVY\
YKEPFWRKKDYCGTMIIDGEEAPVAYTLDDTKPEGNYAAIMGFILAHKARKLARLTKEER\
LKKLCELYAKVLGSLEALEPVHYEEKNWCEEQYSGGCYTTYFPPGILTQYGRVLRQPVDR\
IYFAGTETATHWSGYMEGAVEAGERAAREILHAMGKIPEDEIWQSEPESVDVPAQPITTT\
FLERHLPSVPGLLRLIGLTTIFSATALGFLAHKRGLLVRV",
    "HMGCR": "MLSRLFRMHGLFVASHPWEVIVGTVTLTICMMSMNMFTGNNKICGWNYECPKFEEDVLSS\
DIIILTITRCIAILYIYFQFQNLRQLGSKYILGIAGLFTIFSSFVFSTVVIHFLDKELTG\
LNEALPFFLLLIDLSRASTLAKFALSSNSQDEVRENIARGMAILGPTFTLDALVECLVIG\
VGTMSGVRQLEIMCCFGCMSVLANYFVFMTFFPACVSLVLELSRESREGRPIWQLSHFAR\
VLEEEENKPNPVTQRVKMIMSLGLVLVHAHSRWIADPSPQNSTADTSKVSLGLDENVSKR\
IEPSVSLWQFYLSKMISMDIEQVITLSLALLLAVKYIFFEQTETESTLSLKNPITSPVVT\
QKKVPDNCCRREPMLVRNNQKCDSVEEETGINRERKVEVIKPLVAETDTPNRATFVVGNS\
SLLDTSSVLVTQEPEIELPREPRPNEECLQILGNAEKGAKFLSDAEIIQLVNAKHIPAYK\
LETLMETHERGVSIRRQLLSKKLSEPSSLQYLPYRDYNYSLVMGACCENVIGYMPIPVGV\
AGPLCLDEKEFQVPMATTEGCLVASTNRGCRAIGLGGGASSRVLADGMTRGPVVRLPRAC\
DSAEVKAWLETSEGFAVIKEAFDSTSRFARLQKLHTSIAGRNLYIRFQSRSGDAMGMNMI\
SKGTEKALSKLHEYFPEMQILAVSGNYCTDKKPAAINWIEGRGKSVVCEAVIPAKVVREV\
LKTTTEAMIEVNINKNLVGSAMAGSIGGYNAHAANIVTAIYIACGQDAAQNVGSSNCITL\
MEASGPTNEDLYISCTMPSIEIGTVGGGTNLLPQQACLQMLGVQGACKDNPGENARQLAR\
IVCGTVMAGELSLMAALAAGHLVKSHMIHNRSKINLQDLQGACTKKTA",
    "SULT1A3": "MELIQDTSRPPLEYVKGVPLIKYFAEALGPLQSFQARPDDLLINTYPKSGTTWVSQILDMIYQGGDLEKCNRAPIYVRVPFLEVNDPGEPSGLETLKDTPPPRLIKSHLPLALLPQTLLDQKVKVVYVARNPKDVAVSYYHFHRMEKAHPEPGTWDSFLEKFMAGEVSYGSWYQHVQEWWELSRTHPVLYLFYEDMKENPKREIQKILEFVGRSLPEETMDFMVQHTSFKEMKKNPMTNYTTVPQELMDHSISPFMRKGMAGDWKTTFTVAQNERFDADYAEKMAGCSLSFRSEL",
    "ADRB1": "MGAGVLVLGASEPGNLSSAAPLPDGAATAARLLVPASPPASLLPPASESPEPLSQQWTAG\
MGLLMALIVLLIVAGNVLVIVAIAKTPRLQTLTNLFIMSLASADLVMGLLVVPFGATIVV\
WGRWEYGSFFCELWTSVDVLCVTASIETLCVIALDRYLAITSPFRYQSLLTRARARGLVC\
TVWAISALVSFLPILMHWWRAESDEARRCYNDPKCCDFVTNRAYAIASSVVSFYVPLCIM\
AFVYLRVFREAQKQVKKIDSCERRFLGGPARPPSPSPSPVPAPAPPPGPPRPAAAAATAP\
LANGRAGKRRPSRLVALREQKALKTLGIIMGVFTLCWLPFFLANVVKAFHRELVPDRLFV\
FFNWLGYANSAFNPIIYCRSPDFRKAFQGLLCCARRAARRRHATHGDRPRASGCLARPGP\
PPSPGAASDDDDDDVVGATPPARLLEPWAGCNGGAAADSDSSLDEPCRPGFASESKV",
    "ADRB2": "MGQPGNGSAFLLAPNGSHAPDHDVTQERDEVWVVGMGIVMSLIVLAIVFGNVLVITAIAK\
FERLQTVTNYFITSLACADLVMGLAVVPFGAAHILMKMWTFGNFWCEFWTSIDVLCVTAS\
IETLCVIAVDRYFAITSPFKYQSLLTKNKARVIILMVWIVSGLTSFLPIQMHWYRATHQE\
AINCYANETCCDFFTNQAYAIASSIVSFYVPLVIMVFVYSRVFQEAKRQLQKIDKSEGRF\
HVQNLSQVEQDGRTGHGLRRSSKFCLKEHKALKTLGIIMGTFTLCWLPFFIVNIVHVIQD\
NLIRKEVYILLNWIGYVNSGFNPLIYCRSPDFRIAFQELLCLRRSSLKAYGNGYSSNGNT\
GEQSGYHVEQEKENKLLCEDLPGTEDFVGHQGTVPSDNIDSQGRNCSTNDSLL",
    "DRD2": "MDPLNLSWYDDDLERQNWSRPFNGSDGKADRPHYNYYATLLTLLIAVIVFGNVLVCMAVS\
REKALQTTTNYLIVSLAVADLLVATLVMPWVVYLEVVGEWKFSRIHCDIFVTLDVMMCTA\
SILNLCAISIDRYTAVAMPMLYNTRYSSKRRVTVMISIVWVLSFTISCPLLFGLNNADQN\
ECIIANPAFVVYSSIVSFYVPFIVTLLVYIKIYIVLRRRRKRVNTKRSSRAFRAHLRAPL\
KGNCTHPEDMKLCTVIMKSNGSFPVNRRRVEAARRAQELEMEMLSSTSPPERTRYSPIPP\
SHHQLTLPDPSHHGLHSTPDSPAKPEKNGHAKDHPKIAKIFEIQTMPNGKTRTSLKTMSR\
RKLSQQKEKKATQMLAIVLGVFIICWLPFFITHILNIHCDCNIPPVLYSAFTWLGYVNSA\
VNPIIYTTFNIEFRKAFLKILHC",
}

def get_job_name():
    '''
    read in job name from command line argument or input prompt

    Jobname should be a text file with the job name in it, e.g. "my_job_name.txt"
    structure should be: 
    protein: [sequence in one letter format]
    ligand_smiles: [SMILES string]
    (below are optional, but if provided, they will be used to guide the geometry of the predicted structure)
    contact_residues: [comma separated list of residue indices, 0-indexed]
    max_distance: [maximum distance from the pocket residues to the binder, in angstroms]
    close_residues: [comma separated list of two residue indices, 0-indexed]
    res_distance: [maximum distance between the two residues, in angstroms]
    template_type: [type of template, either 'pdb' or 'url']
    template_pdbid: [PDB ID of the template structure, if template_type is 'url']
    template_file: [path to the template file, if template_type is 'pdb']

    currently only supports one protein chain and one ligand. 
    '''

    job_name = sys.argv[1] if len(sys.argv) > 1 else None
    if job_name is None:
        job_name = input('enter job name: ')

    if not job_name.endswith('.txt'):
        job_name += '.txt'

    print('----------------------------------------------------')
    print(f'Boltz job name: {job_name.replace(".txt", "")}')
    return job_name


def parse_input(job_name: str):
    
    '''read in protein sequence, ligand smiles, contact residues (optional), and 
    max distance (optional) from job name text file, then
    check if protein sequence is an actual sequence or a key in the sequence bank, 
    if it is a key, replace it with the actual sequence.

    for the smiles, check if it is a valid smiles string using regex, if not, 
    query pubchem for the smiles string using the name. If not found, exit with error message.

    read in template type, template pdbid, and template file (optional) from job name text file, then
    check if template type is 'url' or 'pdb', if 'url', check if template pdbid is provided, if 'pdb', 
    check if template file is provided.
    '''

    with open(job_name, 'r') as f:
        lines = f.readlines()
        protein_sequence = lines[0].strip().split(':')[1].strip()
        ligand_smiles = lines[1].strip().split(':')[1].strip()

        contact_residues = []
        max_distance = 5.0
        close_residues = []
        res_distance = 5.0
        template_type = None
        template_pdbid = None
        template_file = None

        if len(lines) > 2:
            for line in lines[2:]:
                if line.startswith('contact_residues'):
                    contact_residues = [int(x) for x in line.strip().split(':')[1].strip().split(',')]
                elif line.startswith('max_distance'):
                    max_distance = float(line.strip().split(':')[1].strip())
                elif line.startswith('close_residues'):
                    close_residues = [int(x) for x in line.strip().split(':')[1].strip().split(',')]
                elif line.startswith('res_distance'):
                    res_distance = float(line.strip().split(':')[1].strip())
                elif line.startswith('template_type'):
                    template_type = line.strip().split(':')[1].strip()
                elif line.startswith('template_pdbid'):
                    template_pdbid = line.strip().split(':')[1].strip()
                elif line.startswith('template_file'):
                    template_file = line.strip().split(':')[1].strip()

    if protein_sequence.split(',')[0] == 'bank' and protein_sequence.split(',')[1] in sequence_bank.keys():
        protein_sequence = sequence_bank[protein_sequence.split(',')[1]]
        print(f'Protein sequence found in sequence bank')

    #check if ligand_smiles is a valid SMILES string using regex using imported smiles pattern from sequence_bank.py
    if not re.match(smiles_pattern, ligand_smiles):
        print(f'getting SMILES string from PubChem for {ligand_smiles}')
        ligand_smiles = smiles_node(ligand_smiles)
        if ligand_smiles == 'unknown':
            print(f'Could not find SMILES string for {ligand_smiles} in PubChem')
            sys.exit(1)

    if contact_residues != []:
        print(f'Contact residues: {contact_residues}')
    print(f'Max distance: {max_distance} angstroms')
    if close_residues != []:
        print(f'Close residues: {close_residues}')
    print(f'Residue distance: {res_distance} angstroms')

    return protein_sequence, ligand_smiles, contact_residues, max_distance, close_residues, res_distance, template_type, template_pdbid, template_file

def make_contraints_template(template_type, template_pdbid, template_file, contact_residues, 
                             max_distance, close_residues, res_distance):
    '''make the constraints and template structure for the Boltz API request based 
    on the input parameters'''

    # Template processing logic
    template_structure = None
    if template_type == 'url' and template_pdbid:
        template_structure = {
            'type': 'url',
            'url': f'https://files.rcsb.org/download/{template_pdbid}.pdb'
        }
    elif template_type == 'pdb' and template_file:
        with open(template_file, 'rb') as f:
            encoded = base64.b64encode(f.read()).decode('utf-8')
            template_structure = {
                'type': 'base64',
                'data': encoded,
                'media_type': 'chemical/x-pdb'
            }

    constraints = []
    if contact_residues:
        constraints.append({
            "type": "pocket", # keep the binder near a set of receptor residues
            "binder_chain_id": "B",
            "contact_residues": { "A": contact_residues }, # receptor residues lining the pocket (0-indexed)
            "max_distance_angstrom": max_distance, # distance from the pocket residues to the binder
            "force": False # bias by default; true = hard-enforce
        })
    if close_residues:
        constraints.append({
            "type": "contact", # keep two residues/atoms within a distance
            "token1": { "type": "polymer_contact", "chain_id": "A", "residue_index": close_residues[0] },
            "token2": { "type": "polymer_contact", "chain_id": "A", "residue_index": close_residues[1] },
            "max_distance_angstrom": res_distance, # distance between the two residues
            "force": False
    })
    return constraints, template_structure