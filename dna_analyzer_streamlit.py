# to run -> streamlit run dna_analyzer_streamlit.py

import streamlit as st 
from streamlit_extras.add_vertical_space import add_vertical_space 
import re
import pandas as pd

st.set_page_config(page_title="DNA/RNA Analyzer")


def main():
    # header
    st.title("DNA/RNA Analyzer", anchor=False)
    st.markdown("Get standard bioinformatics stats, and more, for your DNA or RNA sequence.")

    # tabs
    tab1, tab2 = st.tabs(["Primary", "Translation"])

    with tab1:
        # this is for the "Clear" btn to work
        # See "Option 3" here https://docs.streamlit.io/library/advanced-features/button-behavior-and-examples
        begin = st.container() 
        
        # buttons
        col1, col2 = st.columns([1, 7.5])
        with col1:
            st.button("Submit", type="primary")
        with col2:
            if st.button("Clear", type="secondary"):
                st.session_state.input_seq = ""
        
        # input box
        begin.text_area("Enter DNA/RNA Sequence", key="input_seq")
        input_seq = st.session_state.input_seq
        
        # input validation
        processed_seq = preprocess(input_seq)
        if validate(processed_seq) is False:
            st.markdown(''':red[Invalid DNA/RNA sequence was entered.]''')
            # st.error('This is an error', icon="🚨")
        
        # Initialize output variables
        sequence_length = ""
        gc_content = ""
        reverse_seq = ""
        complement_seq = ""
        reverse_complement_seq = ""
        rna_seq = ""
        dna_seq = ""
        protein_sequence = ""
        protein_length = ""
    
        # Update output variables if the input is validated
        # ie don't show anything in output boxes until we know the input sequence is valid
        if validate(processed_seq):
            sequence_length = count_seq_length(processed_seq)
            gc_content = count_gc(processed_seq)
            reverse_seq = reverse(processed_seq)
            complement_seq = complement(processed_seq)
            reverse_complement_seq = reverse(complement_seq)
            rna_seq = transcribe(processed_seq)
            dna_seq = reverse_transcribe(processed_seq)
            protein_sequence = translate(rna_seq)
            protein_length = count_protein_length(protein_sequence)
    
        # Output boxes
        add_vertical_space(2)
        st.text_input("Sequence Length", value=sequence_length)
        st.text_input("GC Content", value=gc_content)
        st.text_area("Reverse", value=reverse_seq)
        st.text_area("Complement", value=complement_seq)
        st.text_area("Reverse Complement", value=reverse_complement_seq)
        if is_rna(processed_seq) == True:
            st.text_area("DNA Sequence", value=dna_seq)
        elif is_rna(processed_seq) == False:
            st.text_area("RNA Sequence", value=rna_seq)
    
        # credit in a wannabe footer
        st.caption("Built by [David Smehlik](https://www.linkedin.com/in/dvsmehlik/) in 2023 as the final project for "
            "[CS50P](https://pll.harvard.edu/course/cs50s-introduction-programming-python).")

    # tab2: "Translation"
    with tab2: 
        # protein sequence box
        st.text_area("Translated Protein Sequence", value=protein_sequence)
        if "_" in protein_sequence:
            st.markdown('"_" at the end represents 1 or 2 nucleotides that are not forming an amino acid.')
        
        # protein length count
        st.text_input("Protein Length", value=protein_length)
        add_vertical_space(2)
        
        # Plot the amino acid frequency
        ami_aci_sequence = aa_to_ami_aci(protein_sequence)
        if validate(processed_seq):
            st.markdown("**Amino Acid Frequency**")
            st.bar_chart(df_amino_frequency(ami_aci_sequence), x="Amino Acid", y="Frequency", color=(255, 75, 75))




# FUNCTIONS

def preprocess(seq):
    # Remove spaces, new lines, and convert to uppercase
    processed_seq = seq.replace(' ', '').replace('\n', '').upper()
    return processed_seq


def validate(seq):
    if seq:
        # Define a regex pattern for valid DNA characters (A, T, C, G, U)
        pattern = "^[ATCGU]+$"
        
        # Check if the DNA sequence matches the pattern
        if re.match(pattern, seq):
            return True
        else:
            return False


def is_rna(seq):
    if seq:
        if "U" in seq:
            return True
        elif "T" in seq:
            return False


def count_seq_length(seq):
    if seq:
        nucleotide_count = len(seq)
        if nucleotide_count == 1:
            return "1 nucleotide"
        else:
            return f"{nucleotide_count} nucleotides"


def count_protein_length(protein_sequence):
    # don't show "0 amino acids" when user hasn't input a seq
    if protein_sequence:
        # don't count "_" and stop codons, as they're not amino acids
        cleaned_protein_sequence = protein_sequence.replace("_", "").replace("*", "")
        amino_acid_count = len(cleaned_protein_sequence)
        if amino_acid_count == 1:
            return "1 amino acid"
        else:
            return f"{amino_acid_count} amino acids"


def count_gc(seq):
    if seq:
        # Count the occurrences of 'G' and 'C'
        gc = sum(1 for base in seq if base in 'GC')
        
        # Calculate GC content as a percentage
        gc = (gc / len(seq)) * 100
        
        return f"{gc:.2f} %"


def reverse(seq):
    return seq[::-1]


def complement(seq):
    # Create a translation table to replace bases
    complement_table = str.maketrans('ATCGU', 'TAGCA')
    
    # Use translate to replace each base with its complement
    complement_seq = seq.translate(complement_table)
    
    return complement_seq


def transcribe(seq):
    rna = seq.replace('T', 'U')
    return rna


def reverse_transcribe(seq):
    dna = seq.replace('U', 'T')
    return dna


def translate(rna):
    # Define the codon table mapping RNA codons to amino acids
    codon_table = {
        'UUU': 'F', 'UUC': 'F', 
        'UUA': 'L', 'UUG': 'L',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 
        'UAU': 'Y', 'UAC': 'Y',
        'UGU': 'C', 'UGC': 'C',
        'UGG': 'W', 
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L', 
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 
        'CAU': 'H', 'CAC': 'H', 
        'CAA': 'Q', 'CAG': 'Q', 
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I',
        'AUG': 'M',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'AAU': 'N', 'AAC': 'N', 
        'AAA': 'K', 'AAG': 'K',
        'AGU': 'S', 'AGC': 'S', 
        'AGA': 'R', 'AGG': 'R',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V', 
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 
        'GAU': 'D', 'GAC': 'D', 
        'GAA': 'E', 'GAG': 'E', 
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
        # stop codons
        'UGA': '*', 'UAA': '*', 'UAG': '*',
    }
    
    # Split the RNA sequence into codons (triplets) in a list
    codons = [rna[i:i + 3] for i in range(0, len(rna), 3)] 
    
    # Translate each codon to its corresponding amino acid
    protein_sequence = [codon_table.get(codon, '_') for codon in codons]
    protein_sequence = "".join(protein_sequence) # join a list, eg: [G, E, D, A, V, *, T, R, _] into str "GEDAV*TR_"
    
    return protein_sequence 


def aa_to_ami_aci(aa_seq):
    # Define a mapping from single-letter amino acid codes to three-letter codes
    aa_mapping = {
        'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp',
        'C': 'Cys', 'Q': 'Gln', 'E': 'Glu', 'G': 'Gly',
        'H': 'His', 'I': 'Ile', 'L': 'Leu', 'K': 'Lys',
        'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser',
        'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
        '*': '*',
    }
    
    # Map the single-letter amino acid codes to three-letter codes
    ami_aci_sequence = [aa_mapping.get(aa, '_') for aa in aa_seq]
    
    return ami_aci_sequence # a list, eg: [Gly, Asp, Tyr, Leu]


def df_amino_frequency(ami_aci_sequence):
    # Filter out stop codons and unknown amino acids
    filtered_sequence = [amino_acid for amino_acid in ami_aci_sequence if amino_acid not in ['*', '_']]
    
    # Create a DataFrame from the filtered sequence
    df = pd.DataFrame(filtered_sequence, columns=['Amino Acid'])
    
    # Count the frequency of each amino acid
    amino_acid_counts = df['Amino Acid'].value_counts()
    
    # Create a new DataFrame with "Amino Acid" and "Frequency" columns
    result_df = pd.DataFrame({
        'Amino Acid': amino_acid_counts.index,
        'Frequency': amino_acid_counts.values
    })
    
    # Sort by frequency, highest to lowest
    result_df = result_df.sort_values(by='Frequency', ascending=False)
    
    return result_df



if __name__ == "__main__":
    main()