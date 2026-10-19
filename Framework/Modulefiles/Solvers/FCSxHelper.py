import ufl
from mpi4py import MPI
import dolfinx
import meshio
from scipy.spatial import cKDTree
from dolfinx.io import XDMFFile

def tensor_to_voigt(tensor, dim=2):
    if dim == 2:
        # Standard 2D Voigt (Plane Stress/Strain)
        return ufl.as_vector([tensor[0, 0],
                              tensor[1, 1],
                              tensor[0, 1]])
    else:
        # Standard 3D Voigt
        return ufl.as_vector([tensor[0, 0],
                              tensor[1, 1],
                              tensor[2, 2],
                              tensor[1, 2],
                              tensor[0, 2],
                              tensor[0, 1]])

def voigt_to_tensor(v, dim=2):
    if dim == 2:
        # Vector: [sigma_xx, sigma_yy, sigma_xy]
        return ufl.as_tensor([[v[0], v[2]],
                              [v[2], v[1]]])
    else:
        # Vector: [sigma_xx, sigma_yy, sigma_zz, sigma_yz, sigma_xz, sigma_xy]
        return ufl.as_tensor([[v[0], v[5], v[4]],
                              [v[5], v[1], v[3]],
                              [v[4], v[3], v[2]]])

def load_field_to_function(domain, filename, field_name, step_idx=0):
    """
    Loads a specific data field from an XDMF using meshio.TimeSeriesReader
    and injects it into a FEniCSx Function.
    """
    # 1. Create the FunctionSpace and Function
    # Most XDMF node data is stored as Lagrange Degree 1 (linear)
    V = dolfinx.fem.functionspace(domain, ("Lagrange", 1))
    u = dolfinx.fem.Function(V)
    u.name = field_name


    # 2. Use TimeSeriesReader to get the raw NumPy data
    try:
        with meshio.xdmf.TimeSeriesReader(filename) as reader:
            points, cells = reader.read_points_cells()
            # Read metadata for the requested step
            t, point_data, cell_data = reader.read_data(step_idx)
            if field_name not in point_data:
                raise KeyError(f"Field '{field_name}' not found in {filename}. "
                               f"Available: {list(point_data.keys())}")

            raw_data = point_data[field_name]
            raw_data = raw_data.flatten()
            xdmf_coords = points[:, :domain.geometry.dim]
    except Exception as e:
        print(f">>> Error reading {filename} with meshio: {e}")
        return u  # Returns zero-initialized function as fallback

    dof_coords = V.tabulate_dof_coordinates()[:, :domain.geometry.dim]
    tree = cKDTree(xdmf_coords)
    _, xdmf_indices = tree.query(dof_coords)

    # 3. Inject data into the FEniCSx Function
    # NOTE: This assumes the mesh in 'domain' has the same node ordering as 'filename'
    # .real handles potential complex-type builds of PETSc
    if u.x.array.shape == raw_data.shape:
        u.x.array[:] = raw_data[xdmf_indices].real
    else:
        # If there is a mismatch (e.g. ghost cells or MPI), we only fill the local nodes
        # This part is crucial for stability in parallel/larger meshes
        local_size = V.dofmap.index_map.size_local
        u.x.array[:local_size] = raw_data[xdmf_indices[:local_size]].real

    # Update ghost values across MPI processes
    u.x.scatter_forward()

    print(f">>> Successfully loaded '{field_name}' from step {step_idx} (t={t:.2e}s)")
    return u

def ufl_piecewise_linear(
    x,
    x_table,
    y_table
):
    """
    Piecewise-linear interpolation of a UFL expression.

    x:
        UFL scalar expression

    x_table:
        monotonically increasing numpy array

    y_table:
        numpy array of corresponding values

    Returns:
        UFL scalar expression
    """

    expr = ufl.FloatValue(
        float(y_table[-1])
    )

    for i in range(
        len(x_table) - 2,
        -1,
        -1
    ):

        x0 = float(x_table[i])
        x1 = float(x_table[i + 1])

        y0 = float(y_table[i])
        y1 = float(y_table[i + 1])

        slope = (
            (y1 - y0)
            / (x1 - x0)
        )

        yi = (
            y0
            + slope * (x - x0)
        )

        expr = ufl.conditional(
            ufl.lt(x, x1),
            yi,
            expr
        )

    return expr

def load_datastream_to_fenicsx(filename="Datastream.xdmf"):
    """
    Reads the mesh from Datastream.xdmf and converts it
    into a DOLFINx (FEniCSx) Mesh/Domain.
    """
    try:
        # 1. Open the file in read mode
        # We use the COMM_WORLD communicator for parallel compatibility
        with XDMFFile(MPI.COMM_WORLD, filename, "r") as xdmf:

            # 2. Read the mesh (domain)
            # 'mesh' is the default name meshio uses in the XDMF header
            domain = xdmf.read_mesh(name="mesh")

        print(f"Successfully loaded {filename} into FEniCSx domain.")
        print(
            f"Nodes: {domain.geometry.x.shape[0]}, Cells: {domain.topology.index_map(domain.topology.dim).size_global}")

        return domain

    except Exception as e:
        print(f"Error loading mesh to FEniCSx: {e}")
        return None