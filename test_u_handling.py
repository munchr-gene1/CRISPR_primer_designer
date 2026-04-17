
from crispr_primer_designer import CutSite, TemplateSequence, Primer3InputGenerator, Primer3Config

def test_u_matching():
    # gRNA with U
    grna_seq = "GGUCCUCGUGGCCUGGUACA"
    # Template with T (DNA version)
    template_seq = "AGCTGGTCCTCGTGGCCTGGTACAGGG"
    
    cut_site = CutSite(
        target_gene="Test",
        grna_name="Test_g1",
        sequence=grna_seq,
        location="1"
    )
    
    template = TemplateSequence(
        name="Test_Template",
        sequence=template_seq
    )
    
    config = Primer3Config()
    generator = Primer3InputGenerator(config)
    
    matches = generator.find_sequence_in_template(cut_site.sequence, template.sequence)
    
    print(f"gRNA sequence: {cut_site.sequence}")
    print(f"Template sequence: {template.sequence}")
    print(f"Matches found: {len(matches)}")
    
    if len(matches) > 0:
        print("SUCCESS: Match found!")
    else:
        print("FAILURE: No match found (expected failure before fix)")

if __name__ == "__main__":
    test_u_matching()
