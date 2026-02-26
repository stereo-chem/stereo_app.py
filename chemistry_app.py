import streamlit as st
import pubchempy as pcp
import requests
from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem.Draw import rdMolDraw2D
from stmol import showmol
import py3Dmol
import base64

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
        res = requests.get(opsin_url, timeout=5)
        if res.status_code == 200: return res.json()['smiles']
    except: pass
    try:
        pcp_res = pcp.get_compounds(name, 'name')
        if pcp_res: return pcp_res[0].isomeric_smiles
    except: pass
    return None

def render_pro_2d(mol):
    """رسم 2D بنظام SVG لضمان الظهور بوضوح على اللاب توب مع الـ Wedges"""
    mc = Chem.Mol(mol)
    rdDepictor.Compute2DCoords(mc)
    Chem.AssignStereochemistry(mc, force=True, cleanIt=True)
    mc = rdMolDraw2D.PrepareMolForDrawing(mc)
    
    # استخدام SVG لتجنب مشاكل الاختفاء في المتصفح
    drawer = rdMolDraw2D.MolDraw2DSVG(500, 500)
    opts = drawer.drawOptions()
    opts.bondLineWidth = 4.0
    opts.addStereoAnnotation = True
    opts.explicitMethyl = True
    opts.fixedBondLength = 35
    
    drawer.DrawMolecule(mc)
    drawer.FinishDrawing()
    
    svg = drawer.GetDrawingText()
    # تحويل الـ SVG إلى Base64 لعرضه كصورة
    b64 = base64.b64encode(svg.encode('utf-8')).decode("utf-8")
    return f"data:image/svg+xml;base64,{b64}"

# --- 3. المعالجة الرئيسية ---
compound_name = st.text_input("Enter Structure Name:", "1,3-Dimethyl-3-phenylallene")

if st.button("Analyze & Visualize Isomers"):
    with st.spinner('Processing structure...'):
        smiles = get_smiles_smart(compound_name)
        if smiles:
            mol = Chem.MolFromSmiles(smiles)
            allene_p = Chem.MolFromSmarts("C=C=C")
            
            # تحضير الألين للتعامل معه كمركز كيرالي (Axial)
            if mol.HasSubstructMatch(allene_p):
                for match in mol.GetSubstructMatches(allene_p):
                    # تمييز الذرات الطرفية للألين
                    mol.GetAtomWithIdx(match[0]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
                    mol.GetAtomWithIdx(match[2]).SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)

            opts = StereoEnumerationOptions(tryEmbedding=True, onlyUnassigned=False)
            isomers = list(EnumerateStereoisomers(mol, options=opts))
            
            # إذا وجد أيزومر واحد للألين، نقوم بتوليد الصورة المرآتية له يدوياً
            if len(isomers) == 1 and mol.HasSubstructMatch(allene_p):
                iso2 = Chem.Mol(isomers[0])
                for a in iso2.GetAtoms():
                    if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
                        new_tag = Chem.ChiralType.CHI_TETRAHEDRAL_CCW if a.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW else Chem.ChiralType.CHI_TETRAHEDRAL_CW
                        a.SetChiralTag(new_tag)
                isomers.append(iso2)

            st.subheader(f"Found {len(isomers)} Stereoisomer(s)")
            
            st.write("---")
            cols = st.columns(len(isomers))
            
            for i, iso in enumerate(isomers):
                with cols[i]:
                    axial_type = "Ra" if i == 0 else "Sa"
                    st.markdown(f"### Isomer {i+1}: <span style='color: #800000;'>{axial_type}</span>", unsafe_allow_html=True)
                    
                    # 1. عرض الـ 2D (نظام SVG)
                    st.image(render_pro_2d(iso), use_container_width=True)
                    
                    # 2. عرض الـ 3D التفاعلي
                    m3d = Chem.AddHs(iso)
                    AllChem.EmbedMolecule(m3d, AllChem.ETKDG())
                    
                    view = py3Dmol.view(width=400, height=300)
                    view.addModel(Chem.MolToMolBlock(m3d), 'mol')
                    
                    # تلوين ذرات الكربون الطرفية في الألين باللون الأحمر
                    allene_matches = iso.GetSubstructMatches(allene_p)
                    terminal_indices = []
                    if allene_matches:
                        for match in allene_matches:
                            terminal_indices.extend([match[0], match[2]])
                    
                    # تطبيق التلوين والستايل
                    for idx in terminal_indices:
                        view.setStyle({'serial': idx + 1}, {'sphere': {'color': '#FF0000', 'scale': 0.4}, 'stick': {'color': '#FF0000'}})
                    
                    view.setStyle({'not': {'serial': [idx + 1 for idx in terminal_indices]}}, {'stick': {}, 'sphere': {'scale': 0.25}})
                    
                    view.zoomTo()
                    showmol(view)
        else:
            st.error("Compound not found. Please check the name.")
