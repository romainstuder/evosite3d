from modeller import *
from modeller.automodel import *

# Set up environment
env = Environ()
env.io.atom_files_directory = [".", "./atom_files"]

# Read parameters
env.libs.topology.read(file="$(LIB)/top_heav.lib")
env.libs.parameters.read(file="$(LIB)/par.lib")


class MyModel(AutoModel):
    def special_restraints(self, aln):
        """Add custom restraints to maintain scaffold structure"""
        rsr = self.restraints
        at = self.atoms

        # Add secondary structure restraints for helices
        # Assuming residues 1-20 and 25-40 are helices
        rsr.add(secondary_structure.Alpha(self.residue_range("1:A", "20:A")))
        rsr.add(secondary_structure.Alpha(self.residue_range("25:A", "40:A")))

        # Add distance restraints to maintain overall fold
        # Constrain CA-CA distances from template
        for i in range(1, 41, 5):  # Every 5th residue
            for j in range(i + 10, min(i + 30, 41), 5):
                rsr.add(
                    forms.Gaussian(
                        group=physical.xy_distance,
                        feature=features.Distance(at[f"CA:{i}:A"], at[f"CA:{j}:A"]),
                        mean=10.0,
                        stdev=2.0,
                    )
                )

    def select_atoms(self):
        """Select atoms to be optimized"""
        # Optimize all atoms
        return Selection(self)


# Set up automodel
a = MyModel(
    env,
    alnfile="alignment.ali",
    knowns="template_ca",
    sequence="target",
    assess_methods=(assess.DOPE, assess.GA341),
)

# Model settings
a.starting_model = 1
a.ending_model = 5  # Generate 5 models
a.md_level = refine.slow  # Thorough optimization

# Build models
a.make()

# Get best model by DOPE score
ok_models = [x for x in a.outputs if x["failure"] is None]
ok_models.sort(key=lambda x: x["DOPE score"])
best_model = ok_models[0]

print(f"\nBest model: {best_model['name']}")
print(f"DOPE score: {best_model['DOPE score']}")
print(f"GA341 score: {best_model['GA341 score'][0]}")
