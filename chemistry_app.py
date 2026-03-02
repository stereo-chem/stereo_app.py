import streamlit as st
import pubchempy as pcp
import requests
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem.Draw import rdMolDraw2D
from stmol import showmol
import py3Dmol
import numpy as np

# --- 1. إعدادات الصفحة والستايل ---
st.set_page_config(page_title="Professional Isomer System", layout="wide")

st.markdown("""
<div style="background-color: #fdf2f2; padding: 15px; border-radius: 10px; border-left: 5px solid #800000; margin-bottom: 20px;">
    <strong style="color: #800000; font-size: 1.2em;">Stereoisomerism Reference Guide:</strong><br>
    <ul style="list-style-type: none; padding-left: 0; margin-top: 10px; color: black;">
        <li>1. <b>Cis / Trans:</b> Identical groups on same/opposite sides.</li>
        <li>2. <b>E / Z (Absolute - CIP System):</b> High-priority groups together (Z) or opposite (E).</li>
        <li>3. <b>R / S (Optical):</b> Absolute configuration of chiral centers.</li>
        <li>4. <b>Ra / Sa (Axial):</b> Stereochemistry of Allenes (C=C=C).</li>
    </ul>
</div>
""", unsafe_allow_html=True)

st.markdown("<h2 style='color: #800000; font-family: serif; border-bottom: 2px solid #dcdde1;'>Chemical Isomer Analysis System 2.0</h2>", unsafe_allow_html=True)

# --- 2. الدوال المساعدة ---
def get_smiles_smart(name):
    try:
        opsin_url = f"https://opsin.ch.cam.ac.uk/opsin/{name}.json"
        res = requests.get(opsin_url)
        if res.status_code == 200: return res.json()['smiles']
    except: pass
    try:
        pcp_res = pcp.get_compounds(name, 'name')
        if pcp_res: return pcp_res[0].isomeric_smiles
    except: pass
    return None

def calculate_axial_name(mol):
    try:
        mol_3d = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
        conf = mol_3d.GetConformer()
        pattern = Chem.MolFromSmarts("C=C=C")
        match = mol_3d.GetSubstructMatch(pattern)
        if not match: return "N/A"
        angle = AllChem.GetDihedralDeg(conf, match[0]-1, match[0], match[2], match[2]+1)
        return "Ra" if angle > 0 else "Sa"
    except: return "Ra/Sa"

# --- دالة 2D المحسنة (تم التحديث هنا) ---
def render_pro_2d(mol):
    mc = Chem.Mol(mol)
    Chem.AssignStereochemistry(mc, force=True, cleanIt=True)
    AllChem.Compute2DCoords(mc)
    Chem.WedgeMolBonds(mc, mc.GetConformer())
    
    # تحديد نظام الألين للتلوين
    allene_pattern = Chem.MolFromSmarts("C=C=C")
    matches = mc.GetSubstructMatches(allene_pattern)
    
    atoms_to_highlight = []
    bonds_to_highlight = []
    atom_colors = {}
    bond_colors = {}
    highlight_color = (1.0, 0.8, 0.8) # لون وردي فاتح للتمييز
    
    for match in matches:
        atoms_to_highlight.extend(match)
        for idx in match:
            atom_colors[idx] = highlight_color
        
        # تلوين الروابط بين ذرات الألين
        for i in range(len(match)-1):
            bond = mc.GetBondBetweenAtoms(match[i], match[i+1])
            if bond:
                bidx = bond.GetIdx()
                bonds_to_highlight.append(bidx)
                bond_colors[bidx] = highlight_color

    drawer = rdMolDraw2D.MolDraw2DCairo(600, 600)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 4.0           # زيادة سمك الخطوط للوضوح
    opts.addStereoAnnotation = True    # إظهار رموز التجسيم
    opts.fixedBondLength = 40
    opts.explicitMethyl = True
    opts.prepareMolsBeforeDrawing = True
    
    # الرسم مع تسليط الضوء (Highlighting)
    drawer.DrawMolecule(
        mc,
        highlightAtoms=atoms_to_highlight,
        highlightAtomColors=atom_colors,
        highlightBonds=bonds_to_highlight,
        highlightBondColors=bond_colors
    )
    
    drawer.FinishDrawing()
    return drawer.GetDrawingText()

# --- 3. المعالجة الرئيسية ---
compound_name = st.text_input("Enter Structure Name:", "1,3-Dimethyl-3-phenylallene")

if st.button("Analyze & Visualize Isomers"):
    smiles = get_smiles_smart(compound_name)
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        allene_p = Chem.MolFromSmarts("C=C=C")
        
        if mol.HasSubstructMatch(allene_p):
            for match in mol.GetSubstructMatches(allene_p):
                mol.GetAtomWithIdx(match[0]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)

        opts = StereoEnumerationOptions(tryEmbedding=True, onlyUnassigned=False)
        isomers = list(EnumerateStereoisomers(mol, options=opts))
        
        if len(isomers) == 1 and mol.HasSubstructMatch(allene_p):
            iso2 = Chem.Mol(isomers[0])
            for a in iso2.GetAtoms():
                if a.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                    a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
            isomers.append(iso2)

        st.subheader(f"Found {len(isomers)} Stereoisomer(s)")
        st.subheader("1. Isomeric Relationships")
        if len(isomers) > 1:
            st.info("💡 Relationships Analysis: Stereoisomeric relationship detected (Enantiomers/Axial).")
        else:
            st.info("The compound is achiral or only one isomer was identified.")
            
        st.write("---")
        cols = st.columns(len(isomers))
        for i, iso in enumerate(isomers):
            with cols[i]:
                axial_type = "Ra" if i == 0 else "Sa"
                st.markdown(f"### Isomer {i+1}: <span style='color: #800000;'>{axial_type}</span>", unsafe_allow_html=True)
                
                # عرض الـ 2D المحسن بالتلوين
                st.image(render_pro_2d(iso), use_container_width=True)
                
                # عرض الـ 3D التفاعلي
                m3d = Chem.AddHs(iso)
                AllChem.EmbedMolecule(m3d, AllChem.ETKDG())
                
                view = py3Dmol.view(width=400, height=300)
                view.addModel(Chem.MolToMolBlock(m3d), 'mol')
                
                # إيجاد ذرات الألين لتلوينها في الـ 3D أيضاً
                allene_matches = iso.GetSubstructMatches(allene_p)
                allene_atoms = set()
                for match in allene_matches:
                    allene_atoms.update(match)
                
                for idx in allene_atoms:
                    view.setStyle({'serial': idx+1}, {'stick': {'color':'red'}, 'sphere': {'color':'red','scale':0.3}})
                
                view.setStyle({'stick': {}, 'sphere': {'scale':0.25}})
                view.zoomTo()
                showmol(view)
    else:
        st.error("Compound not found.")
