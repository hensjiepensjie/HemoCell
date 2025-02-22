/*
 * Name: Hans Pijpers
 * Student ID: 11001658
 * Description:
 * This code will run a small piece of the capillary system of a computer generated mouse brain.
 */

#include "hemocell.h"
#include <helper/voxelizeDomain.h>
#include "rbcHighOrderModel.h"
#include "pltSimpleModel.h"
#include "cellInfo.h"
#include "fluidInfo.h"
#include "particleInfo.h"
#include "writeCellInfoCSV.h"
#include "preInlet.h"
#include <fenv.h>

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace hemo;


/*Main program*/
int main (int argc, char * argv[]) {

  /*Check for config.xml*/
  if(argc < 2) {
    cout << "Usage: " << argv[0] << " <configuration.xml>" << endl;
    return -1;
  }


  HemoCell hemocell(argv[1], argc, argv);
  Config * cfg = hemocell.cfg;

  /*Define parameters for lattice Boltzmann simulation, read from config.xml */
  hlog << "(Stl preinlet) (Geometry) reading and voxelizing STL file " << (*cfg)["domain"]["geometry"].read<string>() << endl;
  MultiScalarField3D<int> *flagMatrix = 0;
  VoxelizedDomain3D<double> * voxelizedDomain = 0;
  getFlagMatrixFromSTL((*cfg)["domain"]["geometry"].read<string>(),
            (*cfg)["domain"]["fluidEnvelope"].read<int>(),
            (*cfg)["domain"]["refDirN"].read<int>(),
            (*cfg)["domain"]["refDir"].read<int>(),
            voxelizedDomain, flagMatrix,
            (*cfg)["domain"]["blockSize"].read<int>());

  hlog << "(Stl preinlet) (Parameters) setting lbm parameters" << endl;

  param::lbm_base_parameters((*cfg));
  param::printParameters();

  hlog << "(PreInlets) creating preInlet" << endl;
  hemocell.preInlet = new hemo::PreInlet(&hemocell,flagMatrix);

  Box3D slice = flagMatrix->getBoundingBox();

/* Some prints to check the values of the inlet slice 

  hlog << "SLICE x0 = " << slice.x0 << endl;
  hlog << "SLICE x1 = " << slice.x1 << endl;
  hlog << "SLICE y0 = " << slice.y0 << endl;
  hlog << "SLICE y1 = " << slice.y1 << endl;
  hlog << "SLICE z0 = " << slice.z0 << endl;
  hlog << "SLICE z1 = " << slice.z1 << endl;
  
  */

  int max_x = slice.x1;
  int max_y = slice.y1;
  int max_z = slice.z1;

  // You can put the inlet inside the domain, or any other side if needed
  /*Coordinate from y-normal
  slice.z0 = 0;
  slice.z1 = 50;
  slice.y0 = 1;
  slice.y1 = 1;
  slice.x0 = 130;
  slice.x1 = 181;
  */

  // Coordinates from z-normal
  slice.z0 = slice.z1-1.0;
  slice.z1 = slice.z0;
  slice.y0 = 0;
  slice.y1 = 58; // 0.15 * slice.y1
  slice.x0 = 129;
  slice.x1 = 188; // 0.6 * slice.x1
  

  /*
  hlog << "SLICE x0 = " << slice.x0 << endl;
  hlog << "SLICE x1 = " << slice.x1 << endl;
  hlog << "SLICE y0 = " << slice.y0 << endl;
  hlog << "SLICE y1 = " << slice.y1 << endl;
  hlog << "SLICE z0 = " << slice.z0 << endl;
  hlog << "SLICE z1 = " << slice.z1 << endl;

  hlog << "max_x = " << max_x << endl;
  hlog << "max_y = " << max_y << endl;
  hlog << "max_z = " << max_z << endl;
  */
  
  int blockSize = (*cfg)["domain"]["blockSize"].read<int>();
  int envelopeWidth = (*cfg)["domain"]["fluidEnvelope"].read<int>();

  // Direction:: -> define the inflow direction (preInlet is on the X negative side)
  /* y-normal reparallelize
  hemocell.preInlet->preInletFromSlice(Direction::Yneg,slice);
  
  hlog << "(Stl preinlet) (Fluid) Initializing Palabos Fluid Field" << endl;
  MultiBlockManagement3D sparseBlockManagement =
  computeSparseManagement(*plb::reparallelize(*flagMatrix, blockSize, blockSize, blockSize), envelopeWidth);
  hemocell.initializeLattice(sparseBlockManagement);
  */
  
  // Original domain decomposition
  hemocell.preInlet->preInletFromSlice(Direction::Zpos,slice);
  
  hlog << "(Stl preinlet) (Fluid) Initializing Palabos Fluid Field" << endl;
  hemocell.initializeLattice(voxelizedDomain->getMultiBlockManagement());


  /*
  bool sparse = false;
  if (sparse) {
      pcout << "Setting simulation domain mask for sparse decomposition..." << endl;
      MultiScalarField3D<int>* flagMatrix = new MultiScalarField3D<int>(Nx, Ny, Nz);
      setToFunction(*flagMatrix, flagMatrix->getBoundingBox(), FlagMaskDomain3D<unsignedshort>(gfData, 1));

      pcout << "Creating sparse representation ..." << endl;

      //Create sparse representation
      MultiBlockManagement3D sparseBlockManagement =
          computeSparseManagement(*plb::reparallelize(*flagMatrix, blockSize, blockSize, blockSize), envelopeWidth);

      lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(sparseBlockManagement,
          defaultMultiBlockPolicy3D().getBlockCommunicator(),
          defaultMultiBlockPolicy3D().getCombinedStatistics(),
          defaultMultiBlockPolicy3D().getMultiCellAccess<T, DESCRIPTOR>(),
          new BackgroundDynamics(omega));
  }
  else {
      lattice = new MultiBlockLattice3D<T, DESCRIPTOR>(Nx, Ny, Nz, new BackgroundDynamics(omega));
  }
  */

    if (!hemocell.partOfpreInlet) {
      hemocell.lattice->periodicity().toggleAll(false);
    }

    // Setting Preinlet creation
    hemocell.preInlet->initializePreInlet();

    hlog << "(Stl preinlet) (Fluid) Setting up boundaries in Palabos Fluid Field" << endl;
    boundaryFromFlagMatrix(hemocell.lattice,flagMatrix,hemocell.partOfpreInlet);

    hemocell.preInlet->createBoundary();

    /*Disable statistics so it runs faster*/
    hemocell.lattice->toggleInternalStatistics(false);

    /*Equilibrate everything*/
    hemocell.latticeEquilibrium(1.,plb::Array<double, 3>(0.,0.,0.));

    // Driving Force
  hemocell.preInlet->calculateDrivingForce();

  hemocell.lattice->initialize();

  // Adding all the cells
  hemocell.initializeCellfield();

  hemocell.addCellType<RbcHighOrderModel>("RBC", RBC_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("RBC", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());
  hemocell.setInitialMinimumDistanceFromSolid("RBC", 1); //Micrometer! not LU

  hemocell.addCellType<PltSimpleModel>("PLT", ELLIPSOID_FROM_SPHERE);
  hemocell.setMaterialTimeScaleSeparation("PLT", (*cfg)["ibm"]["stepMaterialEvery"].read<int>());

    hemocell.setParticleVelocityUpdateTimeScaleSeparation((*cfg)["ibm"]["stepParticleEvery"].read<int>());

    //hemocell.setRepulsion((*cfg)["domain"]["kRep"].read<double>(), (*cfg)["domain"]["RepCutoff"].read<double>());
    //hemocell.setRepulsionTimeScaleSeperation((*cfg)["ibm"]["stepMaterialEvery"].read<int>());

    vector<int> outputs = {OUTPUT_POSITION,OUTPUT_TRIANGLES,OUTPUT_FORCE,OUTPUT_FORCE_VOLUME,OUTPUT_FORCE_BENDING,OUTPUT_FORCE_LINK,OUTPUT_FORCE_AREA,OUTPUT_FORCE_VISC};
    hemocell.setOutputs("RBC", outputs);
    hemocell.setOutputs("PLT", outputs);

    outputs = {OUTPUT_VELOCITY,OUTPUT_DENSITY,OUTPUT_FORCE,OUTPUT_BOUNDARY,OUTPUT_BOUNDARY,OUTPUT_SHEAR_STRESS,OUTPUT_SHEAR_RATE,OUTPUT_STRAIN_RATE};
    hemocell.setFluidOutputs(outputs);

    /* Printing out the values of the outlet slice to see

    Box3D test = hemocell.lattice->getBoundingBox();

    hlog << "TEST x0 = " << test.x0 << endl;
    hlog << "TEST x1 = " << test.x1 << endl;
    hlog << "TEST y0 = " << test.y0 << endl;
    hlog << "TEST y1 = " << test.y1 << endl;
    hlog << "TEST z0 = " << test.z0 << endl;
    hlog << "TEST z1 = " << test.z1 << endl;

    */

    // For the main simulation domain we have to define outlets, angles from behind inlet view
    if (!hemocell.partOfpreInlet) {
      Box3D bb = hemocell.lattice->getBoundingBox();

      Box3D outlet(290, 297, 170, 194, 88, 128); // right, first branch (0.926*max_x, 0.947*max_x, 0.448*max_y, 0.509*max_y, 0.265*max_z, 0.382*max_z) z-normal 
      //Box3D outlet(300, max_x, 210, 250, 175, 195); // y-normal

      OnLatticeBoundaryCondition3D<T, DESCRIPTOR>* boundary = new BoundaryConditionInstantiator3D
          < T, DESCRIPTOR, WrappedZouHeBoundaryManager3D<T, DESCRIPTOR> >();

      boundary->addPressureBoundary0P(outlet, *hemocell.lattice, boundary::density); // y-normal 0P
      setBoundaryDensity(*hemocell.lattice, outlet, 1.0);
      
      Box3D outlet2(243, 292, 270, 343, 87, 91); // right second branch down (0.78*max_x, 0.93*max_x, 0.71*max_y, 0.90*max_y, 0.26*max_z, 0.27*max_z) z-normal 
      //Box3D outlet2(235, 300, 245, 270, 280, 340); // y-normal

      boundary->addPressureBoundary2N(outlet2, *hemocell.lattice, boundary::density); // y-normal 1P
      setBoundaryDensity(*hemocell.lattice, outlet2, 1.0);
      
      Box3D outlet3(128, 132, 342, 381, 93, 119); // right third branch end (0.41*max_x, 0.42*max_x, 0.90*max_y, max_y, 0.28*max_z, 0.353*max_z) z-normal 
      //Box3D outlet3(100, 117, 220, 240, 355, max_z); // y-normal

      boundary->addPressureBoundary0N(outlet3, *hemocell.lattice, boundary::density); // y-normal 0N
      setBoundaryDensity(*hemocell.lattice, outlet3, 1.0);
      
      Box3D outlet4(13, 14, 93, 127, 130, 163); // left lower branch (0.043*max_x, 0.044*max_x, 0.246*max_y, 0.333*max_y, 0.39*max_z, 0.486*max_z) z-normal 
      //Box3D outlet4(0, 15, 170, 201, 94, 128); // y-normal
 
      boundary->addPressureBoundary0N(outlet4, *hemocell.lattice, boundary::density); // y-normal 0N
      setBoundaryDensity(*hemocell.lattice, outlet4, 1.0);

      Box3D outlet5(36, 37, 10, 46, 50, 91); // left upper branch (0.116*max_x, 0.117*max_x, 0.028*max_y, 0.12*max_y, 0.15*max_z, 0.27*max_z) z-normal 
      //Box3D outlet5(7, 34, 245, 295, 10, 45); // y-normal

      boundary->addPressureBoundary0N(outlet5, *hemocell.lattice, boundary::density); // y-normal 0N
      setBoundaryDensity(*hemocell.lattice, outlet5, 1.0);
    }

    //loading the cellfield
    if (not cfg->checkpointed) {
      hemocell.loadParticles();
      hemocell.writeOutput();
    } else {
      hemocell.loadCheckPoint();
    }

    //Restructure atomic blocks on processors when possible
    //hemocell.doRestructure(false); // cause errors(?)

    if (hemocell.iter == 0) {
      pcout << "(Stl preinlet) fresh start: warming up cell-free fluid domain for "  << (*cfg)["parameters"]["warmup"].read<plint>() << " iterations..." << endl;
      for (plint itrt = 0; itrt < (*cfg)["parameters"]["warmup"].read<plint>(); ++itrt) {
          hemocell.lattice->collideAndStream();
      }
    }

    unsigned int tmax = (*cfg)["sim"]["tmax"].read<unsigned int>();
    unsigned int tmeas = (*cfg)["sim"]["tmeas"].read<unsigned int>();
    unsigned int tcheckpoint = (*cfg)["sim"]["tcheckpoint"].read<unsigned int>();
    unsigned int tbalance = (*cfg)["sim"]["tbalance"].read<unsigned int>();
    unsigned int tcsv = (*cfg)["sim"]["tcsv"].read<unsigned int>();


    pcout << "(Stl preinlet) Starting simulation..." << endl;

    while (hemocell.iter < tmax ) {
      //preinlet.update();
      hemocell.iterate();

      if (hemocell.partOfpreInlet) {
        //Set driving force as required after each iteration
        hemocell.preInlet->setDrivingForce();
      }

      hemocell.preInlet->applyPreInlet();

      // Load-balancing! Only enable if PARMETIS build is available
      /*
      if (hemocell.iter % tbalance == 0) {
        if(hemocell.calculateFractionalLoadImbalance() > (*cfg)["parameters"]["maxFlin"].read<double>()) {
          hemocell.doLoadBalance();
          hemocell.doRestructure();
        }
      }
      */
      if (hemocell.iter % tmeas == 0) {
          pcout << "(main) Stats. @ " <<  hemocell.iter << " (" << hemocell.iter * param::dt << " s):" << endl;
          pcout << "\t # of cells: " << CellInformationFunctionals::getTotalNumberOfCells(&hemocell);
          pcout << " | # of RBC: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "RBC");
          pcout << ", PLT: " << CellInformationFunctionals::getNumberOfCellsFromType(&hemocell, "PLT") << endl;
          FluidStatistics finfo = FluidInfo::calculateVelocityStatistics(&hemocell); double toMpS = param::dx / param::dt;
          pcout << "\t Velocity  -  max.: " << finfo.max * toMpS << " m/s, mean: " << finfo.avg * toMpS<< " m/s, rel. app. viscosity: " << (param::u_lbm_max*0.5) / finfo.avg << endl;
          ParticleStatistics pinfo = ParticleInfo::calculateForceStatistics(&hemocell); double topN = param::df * 1.0e12;
          pcout << "\t Force  -  min.: " << pinfo.min * topN << " pN, max.: " << pinfo.max * topN << " pN (" << pinfo.max << " lf), mean: " << pinfo.avg * topN << " pN" << endl;

        // Additional useful stats, if needed
         // finfo = FluidInfo::calculateForceStatistics(&hemocell);
          //Set force as required after this function;
          //setExternalVector(*hemocell.lattice, hemocell.lattice->getBoundingBox(),
                    // DESCRIPTOR<T>::ExternalField::forceBeginsAt,
                   //  hemo::Array<T, DESCRIPTOR<T>::d>(poiseuilleForce, 0.0, 0.0));
         // pcout << "Fluid force, Minimum: " << finfo.min << " Maximum: " << finfo.max << " Average: " << finfo.avg << endl;
          //ParticleStatistics pinfo = ParticleInfo::calculateVelocityStatistics(&hemocell);
          //pcout << "Particle velocity, Minimum: " << pinfo.min << " Maximum: " << pinfo.max << " Average: " << pinfo.avg << endl;

          hemocell.writeOutput();
      }
      
      if (hemocell.iter % tcsv == 0) {
          pcout << "Saving simple mean cell values to CSV at timestep " << hemocell.iter << endl;
          writeCellInfo_CSV(hemocell);
      }

      if (hemocell.iter % tcheckpoint == 0) {
          hemocell.saveCheckPoint();
      }
    }

  pcout << "(main) Simulation finished :) " << endl;

  return 0;

}
