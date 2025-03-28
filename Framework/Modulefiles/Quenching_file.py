from .ModuleStructure_file_new import CalcModule
import meshio

class Quenchingmodule(CalcModule):
    def __init__(self,infile):
        infile = "Inputs/" + infile + ".json"
        super().__init__("Quenching", infile)

    def run(self):
        outstr = ["\n---------------------------------------------------------------------\n",
                  "Quenching module: " + self.inputfile + "\n",
                  "Material type: " + self.minput["Materialtype"],
                  "Phase mixing model: " + self.minput["MixtureModel"],
                  "Quenching type: " + self.minput["quenchmedium"],
                  "Starting temp: ",
                  "Quench temp: " + str(self.minput["quenchtemp"]-273.15),
                  "Quenching time: " + str(self.minput["quenchtime"]) + " seconds",
                  "Phases: ",
                  "Element quad: " + str(2),
                  "\n---------------------------------------------------------------------\n\n"]

        for line in outstr:
            self.writeoutput(line)
            print(line)



        from Framework.Modulefiles.Solvers.ComsolSolver import runComsol

        if not self.check_runcondition():
            print("Using precalculated " + str(self.module) + " simulation")

            with meshio.xdmf.TimeSeriesReader("Datastream_Cache.xdmf") as reader:
                points, cells = reader.read_points_cells()
                pd_list, cd_list, t_list = list(), list(), list()
                for k in range(reader.num_steps):
                    t, point_data, cell_data = reader.read_data(k)
                    t_list.append(t)
                    pd_list.append(point_data)
                    cd_list.append(cell_data)
            with meshio.xdmf.TimeSeriesWriter("Datastream.xdmf") as writer:
                writer.write_points_cells(points, cells)
                for i in range(len(t_list)):
                    writer.write_data(t=t_list[i], point_data=pd_list[i], cell_data=cd_list[i])

        else:
            if self.program == "Comsol":
                print('\nQuenching module')

                if self.program == 'FCSx':
                    print('FeniCSx not implemented')

                    # rundocker(self)
                elif self.program == 'Comsol':
                    print('Using COMSOL for FEM calculation')
                    runComsol(self)
            else:
                raise KeyError(str(self.program) + " not implemented in " + str(self.module) + " module")
