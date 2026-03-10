import RNA


def rnafold_structure(sequence):
    """
    Predict RNA secondary structure using RNAfold.
    """
    fc = RNA.fold_compound(sequence)
    structure, mfe = fc.mfe()

    return {
        "sequence": sequence,
        "structure": structure,
        "mfe": mfe,
    }


def rnacofold_structure(srna, mrna):
    """
    Predict duplex structure using RNAcofold.
    """
    combined = srna + "&" + mrna

    fc = RNA.fold_compound(combined)
    structure, mfe = fc.mfe()

    return {
        "sequence": combined,
        "structure": structure,
        "mfe": mfe,
    }
