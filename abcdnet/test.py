import ROOT
import numpy as np

# Number of events to generate for signal and background
n_events = 2500000
# n_events = 25000

# Define the mean and covariance matrices for the 3D gaussian
mean_bkg = np.array([0, 0, 0])
covariance_bkg = np.array([
    [1, -0.8, 0],
    [-0.8, 1, 0],
    [0, 0, 1]
])*1.5*1.5

mean_sig = np.array([2.5, 2.5, 2])
covariance_sig = np.array([
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, 1]
])*1.5*1.5

# Generate random numbers following a multivariate gaussian distribution
rng = np.random.default_rng()
multigauss_3d_data_bkg = rng.multivariate_normal(mean_bkg, covariance_bkg, n_events)
multigauss_3d_data_sig = rng.multivariate_normal(mean_sig, covariance_sig, n_events)

# Open a ROOT file to store the data
f = ROOT.TFile("multigaussian_3d_bkg.root", "RECREATE")

# Create a TTree and set the branch addresses for the 3 variables
tree = ROOT.TTree("tree", "Multivariate Gaussian 3D Distribution")
x = np.zeros(1, dtype=float)
y = np.zeros(1, dtype=float)
z = np.zeros(1, dtype=float)
tree.Branch('x0', x, 'x/D')
tree.Branch('x1', y, 'y/D')
tree.Branch('x2', z, 'z/D')

# Fill the TTree with the multivariate gaussian data
for i in range(n_events):
    x[0] = multigauss_3d_data_bkg[i][0]
    y[0] = multigauss_3d_data_bkg[i][1]
    z[0] = multigauss_3d_data_bkg[i][2]
    tree.Fill()

# Write the TTree and close the ROOT file
f.Write()
f.Close()

f = ROOT.TFile("multigaussian_3d_sig.root", "RECREATE")

# Create a TTree and set the branch addresses for the 3 variables
tree = ROOT.TTree("tree", "Multivariate Gaussian 3D Distribution")
x = np.zeros(1, dtype=float)
y = np.zeros(1, dtype=float)
z = np.zeros(1, dtype=float)
tree.Branch('x0', x, 'x/D')
tree.Branch('x1', y, 'y/D')
tree.Branch('x2', z, 'z/D')

# Fill the TTree with the multivariate gaussian data
for i in range(n_events):
    x[0] = multigauss_3d_data_sig[i][0]
    y[0] = multigauss_3d_data_sig[i][1]
    z[0] = multigauss_3d_data_sig[i][2]
    tree.Fill()

# Write the TTree and close the ROOT file
f.Write()
f.Close()


print("Multivariate 3D Gaussian data saved to multigaussian_3d.root.")

