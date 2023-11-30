# to run -> streamlit run dna_analyzer_streamlit.py

import streamlit as st 
from streamlit_extras.add_vertical_space import add_vertical_space 
import re
import pandas as pd
import altair as alt

st.set_page_config(page_title="DNA/RNA Analyzer")


def main():
    # header
    st.title("DNA/RNA Analyzer", anchor=False)
    st.markdown("Get standard bioinformatics stats, and more, for your DNA or RNA sequence.")
    
    ### TABS ###
    tab1, tab2 = st.tabs(["Primary", "Translation"])
    
    with tab1:
        # this is for the "Clear" btn to work
        # See "Option 3" here https://docs.streamlit.io/library/advanced-features/button-behavior-and-examples
        begin = st.container() 
        
        # Buttons
        col1, col2 = st.columns([0.3, 0.7])
        with col1:
            col1, col2 = st.columns([0.5, 0.5])
            with col1:
                st.button("Submit", type="primary", use_container_width=True)
            with col2:
                if st.button("Clear", type="secondary", use_container_width=True):
                    st.session_state.input_seq = ""
        
        
        # input box
        begin.text_area("Enter DNA/RNA Sequence", key="input_seq")
        input_seq = st.session_state.input_seq
        
        # input validation
        processed_seq = preprocess(input_seq)
        if validate(processed_seq) is False:
            st.markdown(":rotating_light: :red[Invalid characters in entered sequence. Only A, C, G, T, U, including lowercase, "
                "and spaces are allowed.]")
            # st.error('This is an error', icon="🚨")
        
        # show warning msg if both T and U in input
        if "T" in processed_seq and "U" in processed_seq:
            st.markdown(":warning: Entered sequence contains both T and U.")
        
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
        if is_rna(processed_seq) is True:
            st.text_area("DNA Sequence", value=dna_seq)
        elif is_rna(processed_seq) is False:
            st.text_area("RNA Sequence", value=rna_seq)
        
        
        # wide mode msg
        # st.caption("For the app to take the entire width of the screen, go to: Top right → 3 dots menu → Settings → Wide mode")
        # st.markdown("For the app to take the entire width of the screen, go to: Top right → 3 dots menu → Settings → Wide mode")
        
        # credit in a wannabe footer + wide mode msg
        st.caption(":gray[For the app to take the entire width of the screen, go to: Top right → 3 dots menu → Settings → Wide mode.]  \n"
                "Built by [David S.](https://www.linkedin.com/in/davess/) in 2023, originally as the final project for "
                "[CS50P](https://pll.harvard.edu/course/cs50s-introduction-programming-python).")
    
    
    ## Tab2 "Translation" ##
    with tab2: 
        # protein sequence box
        st.text_area("Translated Protein Sequence", value=protein_sequence)
        if "_" in protein_sequence:
            st.markdown('"_" at the end represents 1 or 2 nucleotides that are not forming a codon.')
            # add_vertical_space(1)
        
        
        # codon count
        codon_count = ""
        if protein_sequence:
            codon_count = count_codons(get_codons(rna_seq))
        st.text_input("Codon Count", value=codon_count)
        # add_vertical_space(1)
        
        # show Codon Usage table
        st.markdown("**Codon Usage**")
        st.dataframe(
            create_df_codons(rna_seq), 
            # the following _config is for showing the progress chart in each row
            column_config={
                "Relative Usage": st.column_config.ProgressColumn(
                    # "Relative Usage",
                    # help="The sales volume in USD",
                    # format="%.2f",
                    min_value=0,
                    max_value=1,
                ),
            },
            use_container_width=True,
            hide_index=True,
        )
        # add_vertical_space(1)
        
        
        # "Amino Acid Count" box
        st.text_input("Amino Acid Count", value=protein_length)
        # add_vertical_space(1)
        
        if protein_sequence and protein_length != "0 amino acids":
            # Plot the amino acid frequency
            ami_aci_sequence = aa_to_ami_aci(protein_sequence)      
            result_df = df_amino_frequency(ami_aci_sequence)
            chart = plot_amino_frequency(result_df)
            
            st.markdown("**Amino Acid Frequency**")
            # st.bar_chart(df_amino_frequency(ami_aci_sequence), x="Amino Acid", y="Frequency", color=(255, 75, 75))
            st.altair_chart(chart, use_container_width=True, theme="streamlit")



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


def plot_amino_frequency(df):
    chart = (
        alt.Chart(df)
        .mark_bar()
        .encode(
            x=alt.X('Amino Acid', sort=None),
            y=alt.Y('Frequency', axis=alt.Axis(format='d')),
            size=alt.value(25),
            color=alt.value('rgb(255, 75, 75)'),
            tooltip=['Amino Acid', 'Frequency']
        )
        .configure_legend(disable=True)  # Hide the legend
    )
    
    return chart


def get_codons(rna):
    ## Remove a possible extra 1 or 2 nucleotides, so that they don't show up the codon_list ##
    # Determine the length of the RNA sequence
    seq_length = len(rna)
    
    # Calculate the number of complete codons (triplets)
    num_complete_codons = seq_length // 3
    
    # Remove the extra letters at the end if the length is not divisible by 3
    rna = rna[:num_complete_codons * 3]
    
    # Split the RNA sequence into codons (triplets) in a list
    codon_list = [rna[i:i + 3] for i in range(0, len(rna), 3)] 
    
    return codon_list


def count_codons(codon_list):
    # don't show "0 codons" when user hasn't input a seq
    # if codon_list:
    codon_count = len(codon_list)
    if codon_count == 1:
        return "1 codon"
    else:
        return f"{codon_count} codons"


def codon_to_ami_aci(codon_list):
    # Define the codon table mapping RNA codons to ami aci
    codon_table = {
        'UUU': 'Phe', 'UUC': 'Phe', 
        'UUA': 'Leu', 'UUG': 'Leu',
        'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser', 
        'UAU': 'Tyr', 'UAC': 'Tyr',
        'UGU': 'Cys', 'UGC': 'Cys',
        'UGG': 'Trp', 
        'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu', 
        'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro', 
        'CAU': 'His', 'CAC': 'His', 
        'CAA': 'Gln', 'CAG': 'Gln', 
        'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg', 
        'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile',
        'AUG': 'Met',
        'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
        'AAU': 'Asn', 'AAC': 'Asn', 
        'AAA': 'Lys', 'AAG': 'Lys',
        'AGU': 'Ser', 'AGC': 'Ser', 
        'AGA': 'Arg', 'AGG': 'Arg',
        'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val', 
        'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala', 
        'GAU': 'Asp', 'GAC': 'Asp', 
        'GAA': 'Glu', 'GAG': 'Glu', 
        'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly',
        # stop codons
        'UGA': '*', 'UAA': '*', 'UAG': '*',
    }
    
    # # Create a list of dictionaries for each codon and its corresponding amino acid
    # data = [{'Codon': codon, 'Amino Acid': codon_table.get(codon, 'Unknown')} for codon in codon_list]
    
    
    # Initialize a dictionary to store codon counts
    codon_counts = {}
    
    # Iterate through the codon list and count occurrences
    for codon in codon_list:
        codon_counts[codon] = codon_counts.get(codon, 0) + 1
    
    # Create a list of dictionaries for each codon and its corresponding amino acid and count
    data = [{
        'Codon': codon, 
        'Amino Acid': codon_table.get(codon, 'Unknown'), 
        'Total Count': count,
    } for codon, count in codon_counts.items()]
    
    # Add codons from codon_table that are not already in the data with a count of 0
    for codon, amino_acid in codon_table.items():
        if codon not in codon_counts:
            data.append({
                'Codon': codon,
                'Amino Acid': amino_acid,
                'Total Count': 0,
            })
    
    return data


def create_df_codons(rna):
    # Get the codon list
    codon_list = get_codons(rna)
    
    # Get the data for the DataFrame
    data = codon_to_ami_aci(codon_list)
    
    # Create a DataFrame
    df = pd.DataFrame(data, columns=['Amino Acid', 'Codon', 'Total Count'])
    
    # # Group by codon and count occurrences
    # grouped_df = df.groupby(['Amino Acid', 'Codon']).size().reset_index(name='Count')
    
    
    # Get the total count of each amino acid
    total_count_for_aa = df.groupby('Amino Acid')['Total Count'].transform('sum')
    
    # Calculate relative usage (percentage) for each codon, fill NaN with 0
    df['Relative Usage'] = (df['Total Count'] / total_count_for_aa).fillna(0)
    
    # # Add "%" to the "Relative Usage" column
    # df['Relative Usage'] = df['Relative Usage'].astype(str) + '%'
    # doesn't work bcs streamlit must have only number to create the progress bar chart
    
    
    # Calculate total usage for each codon
    total_count_all_codons = df['Total Count'].sum()
    # df['Total Usage'] = ((df['Total Count'] / total_count_all_codons) * 100).fillna(0).round(2)
    df['Total Usage'] = ((df['Total Count'] / total_count_all_codons) * 100).fillna(0).apply(lambda x: '{:.2f}'.format(x))
    
    # Add "%" to the "Total Usage" column
    df['Total Usage'] = df['Total Usage'].astype(str) + '%'
    
    
    # Rearrange columns
    df = df.reindex(['Amino Acid', 'Codon', 'Total Count', 'Total Usage', 'Relative Usage'], axis=1)
    
    # Sort by 'Amino Acid' and then reverse sort by 'Relative Usage'
    df = df.sort_values(by=['Amino Acid', 'Relative Usage'], ascending=[True, False])
    
    # Apply styling to align the "Total Count" column to the left
    # df = df.style.set_properties(**{'Total Count': 'text-align: left'})
    df['Total Count'] = df['Total Count'].astype(str)
    
    return df



if __name__ == "__main__":
    main()
