import sys
from PyQt6.QtWidgets import (
    QApplication,
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QPushButton,
    QVBoxLayout,
    QWidget,
)
import pyvista as pv
from pyvistaqt import QtInteractor


class XDMFViewerWindow(QMainWindow):

  def __init__(self):
    super().__init__()
    self.setWindowTitle("XDMF Results Viewer")
    self.resize(1100, 750)

    self.dataset = None
    self.current_mesh = None

    # Main layout setup
    central_widget = QWidget(self)
    self.setCentralWidget(central_widget)
    main_layout = QHBoxLayout(central_widget)

    # Left Control Panel
    control_layout = QVBoxLayout()

    # 1. File Load Button
    self.load_btn = QPushButton("Open XDMF File", self)
    self.load_btn.clicked.connect(self.load_file)
    control_layout.addWidget(self.load_btn)

    # 2. Field Selector (Scalars)
    control_layout.addWidget(QLabel("Select Data Field:"))
    self.field_combo = QComboBox(self)
    self.field_combo.currentIndexChanged.connect(self.update_plot)
    control_layout.addWidget(self.field_combo)

    # 3. Block/Time Step Selector
    control_layout.addWidget(QLabel("Select Time Step / Block:"))
    self.block_combo = QComboBox(self)
    self.block_combo.currentIndexChanged.connect(self.update_plot)
    control_layout.addWidget(self.block_combo)

    control_layout.addStretch()

    # Right 3D Plotter Widget (PyVista Qt Interactor)
    self.plotter = QtInteractor(self)
    self.plotter.add_axes()
    self.plotter.show_grid()

    # Add layouts to main window
    main_layout.addLayout(control_layout, stretch=1)
    main_layout.addWidget(self.plotter.interactor, stretch=4)

  def load_file(self):
    """Opens a file dialog to select and load an XDMF file."""
    file_path, _ = QFileDialog.getOpenFileName(
        self, "Open XDMF File", "", "XDMF Files (*.xdmf *.xmf)"
    )
    if not file_path:
      return

    try:
      self.dataset = pv.read(file_path)
      print(f"Loaded: {file_path}")
      self.populate_controls()
      self.update_plot()
    except Exception as e:
      print(f"Error loading file: {e}")

  def populate_controls(self):
    """Populates dropdown menus based on the loaded dataset."""
    self.field_combo.blockSignals(True)
    self.block_combo.blockSignals(True)

    self.field_combo.clear()
    self.block_combo.clear()

    # Extract available fields (scalars/vectors)
    sample_mesh = self.get_active_mesh()
    if sample_mesh and hasattr(sample_mesh, "point_data"):
      fields = list(sample_mesh.point_data.keys()) + list(
          sample_mesh.cell_data.keys()
      )
      self.field_combo.addItems(fields)

    # Handle MultiBlock or Time series dataset structure
    if isinstance(self.dataset, pv.MultiBlock):
      for i in range(len(self.dataset)):
        name = self.dataset.get_block_name(i) or f"Block {i}"
        self.block_combo.addItem(name, userData=i)
    else:
      self.block_combo.addItem("Default Dataset", userData=0)

    self.field_combo.blockSignals(False)
    self.block_combo.blockSignals(False)

  def get_active_mesh(self):
    """Retrieves the currently selected sub-mesh or block."""
    if self.dataset is None:
      return None

    if isinstance(self.dataset, pv.MultiBlock):
      idx = self.block_combo.currentData()
      if idx is not None and idx < len(self.dataset):
        return self.dataset[idx]
      return None
    else:
      return self.dataset

  def update_plot(self):
    """Clears and re-renders the 3D canvas with the selected field and block."""
    if self.dataset is None:
      return

    mesh = self.get_active_mesh()
    if mesh is None or mesh.n_points == 0:
      return

    field_name = self.field_combo.currentText()

    # Clear previous actors from the scene
    self.plotter.clear()
    self.plotter.add_axes()
    self.plotter.show_grid()

    # Add mesh with selected field mapping
    if field_name:
      self.plotter.add_mesh(
          mesh, scalars=field_name, show_edges=True, cmap="jet"
      )
    else:
      self.plotter.add_mesh(mesh, show_edges=True)

    self.plotter.reset_camera()
    self.plotter.render()


if __name__ == "__main__":
  app = QApplication(sys.argv)
  viewer = XDMFViewerWindow()
  viewer.show()
  sys.exit(app.exec())