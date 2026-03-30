# Mutation Annotation Guidelines

The example (phenotypes_example.txt) combines elements from both the *E. coli* and *Salmonella* mutation datasets. This should provide representative examples of the different mutation types we work with.

The annotations are just flat text files in a tab seperated format. Some prefer to open the file in spreadsheet editors like Excel in order to work with them. That is fine, just remember to save the files in the correct format (Tab-delimited Text).

Most of the columns should be self-explanatory, but please reach out if anything requires clarification. Here we mainly explain the **final column**, which is specifically intended for phenotypes that only manifest in combination with other mutations.

For example:

- If mutation **A** produces phenotype **X**, but mutation **A + B** produces phenotype **Y**, these must be recorded as **two separate entries** in the sheet.

---

## Final Column: Format

The final column specifies the additional mutations required for the phenotype to occur. These are annotated using the following format:
```
<GENE NAME><VARIANT NO><ACCESSION NO>_<WT AA><POSITION><MUTANT AA>.<ALTERNATIVE AA>.<ALTERNATIVE AA>
```

### Example
```
gyrA_1_MH933946.1_S83Y.F.A
```

### Explanation
- **First section**
`<GENE NAME><VARIANT NO><ACCESSION NO>`: Must exactly match one of the existing **Gene accession** entries in the first column.
	- Example: `gyrA_1_MH933946.1`

- **Second section**
`<WT AA><POSITION><MUTANT AA>.<ALTERNATIVE AA>.<ALTERNATIVE AA>`: Represents the amino acid substitution(s).
	- Example: `S83Y.F.A`
	- The letters separated by periods indicate **possible alternatives** to the wild-type residue.
	- This may consist of a **single amino acid** or **multiple alternatives**, as shown above.

---

## Multiple Mutations

If a phenotype requires **more than one additional mutation**, separate them with commas.

**Example:**
```
pmrA_1_CP047010.1_S39I,pmrA_1_CP047010.1_R81S
```

In this case, the phenotype occurs only if substitutions are present at **both positions 39 and 81**, in addition to the specified mutation.

---

## Alternative Combinations

In cases where a phenotype may depend on **different combinations of mutations**, separate the options with semicolons.

**Example:**
```
gyrA_1_MH933946.1_S83Y.F.A;gyrA_1_MH933946.1_D87N.G.Y.K
```

Here, the phenotype is observed if, in addition to the specified mutation, a substitution is present at **either position 83 or position 87**.

---

## Duplication of Dependent Mutations

Mutations that depend on one another will typically appear **multiple times** in the dataset.
- Example: There will be one entry for **position 83 dependent on position 87**, and another entry for **position 87 dependent on position 83**.

---

## Notes

If you encounter additional mechanisms or dependencies not accounted for in this framework, please bring them to our attention so they can be incorporated.

## Known issues
- No way to indicate that any disruption of a gene in question will cause resistance (Klebsiella example).
