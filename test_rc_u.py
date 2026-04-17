
from crispr_primer_designer import Primer3InputGenerator, Primer3Config

def test_rc_u():
    config = Primer3Config()
    generator = Primer3InputGenerator(config)
    
    # RNA: GGU (GGU)
    # DNA: GGT (GGT)
    # RC of GGT: ACC
    
    seq_u = "GGU"
    rc = generator._reverse_complement(seq_u)
    print(f"Sequence: {seq_u}")
    print(f"Reverse complement: {rc}")
    
    if rc == "ACC":
        print("SUCCESS: RC of GGU is ACC")
    else:
        print(f"FAILURE: RC of GGU expected ACC but got {rc}")

if __name__ == "__main__":
    test_rc_u()
