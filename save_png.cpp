#include <fstream>
#include <iomanip>
#include <limits>
#include "meshvectorfield.hpp"
#include "dxf_solid.hpp"
#include "mydxffile.hpp"
#include "gtkplotter.hpp"
#include "geomplotter.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "trajectorydiagnostics.hpp"
#include "particledatabase.hpp"
#include "particlediagplotter.hpp"
#include <sstream>



using namespace std;



void simu( int argc, char **argv )
{
    double y_mm = atof( argv[1] ) + 0.0;
    double V    = atof( argv[2] );

    ostringstream tagstream;
    tagstream << y_mm << "_" << V;
    string tag = tagstream.str();

    string geomfile = "geom_" + tag + ".dat";
    string epotfile = "epot_" + tag + ".dat";
    string pdbfile  = "pdb_"  + tag + ".dat";
    string pngfile  = "zy_view_" + tag + ".png";

    std::ifstream is_geom( geomfile.c_str() );
    if( !is_geom.good() )
        throw( Error( ERROR_LOCATION, "couldn\'t open file \'" + geomfile + "\'" ) );
    Geometry geom( is_geom );
    is_geom.close();
    geom.build_surface();

    std::ifstream is_epot( epotfile.c_str() );
    if( !is_epot.good() )
        throw( Error( ERROR_LOCATION, "couldn\'t open file \'" + epotfile + "\'" ) );
    EpotField epot( is_epot, geom );
    is_epot.close();

    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
                                     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
                                     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    std::ifstream is_pdb( pdbfile.c_str() );
    if( !is_pdb.good() )
        throw( Error( ERROR_LOCATION, "couldn\'t open file \'" + pdbfile + "\'" ) );
    ParticleDataBase3D pdb( is_pdb, geom );
    is_pdb.close();

    VectorField *bfield = NULL;
    if( argc >= 4 ) {                               // optional b-field file is now argv[3]
        bool fout[3] = {true, true, true};
        MeshVectorField *mesh_bfield = new MeshVectorField( MODE_3D, fout, 1.0e-3, 1.0, argv[3] );
        field_extrpl_e bfldextrpl[6] = { FIELD_ZERO, FIELD_ZERO,
                                         FIELD_ZERO, FIELD_ZERO,
                                         FIELD_ZERO, FIELD_ZERO };
        mesh_bfield->set_extrapolation( bfldextrpl );
        bfield = mesh_bfield;
    }

    GeomPlotter gplotter( geom );
    gplotter.set_size( 1600, 900 );
    gplotter.set_epot( &epot );
    gplotter.set_particle_database( &pdb );
    if( bfield )
        gplotter.set_bfield( bfield );
    gplotter.set_view_si( VIEW_ZY, 0.0 );
    gplotter.set_ranges( 2.46, -0.01, 2.562, 0.01 );
    gplotter.plot_png( pngfile.c_str() );
}


int main( int argc, char **argv )
{
    if( argc <= 2 ) {
        cerr << "Usage: plot.run <scanner_y_mm> <Vplate_V> <bfield>\n";
        exit( 1 );
    }

    try {
	ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
	ibsimu.set_thread_count( 4 );
	simu( argc, argv );
    } catch( Error e ) {
	e.print_error_message( ibsimu.message( 0 ) );
        exit( 1 );
    }

    return( 0 );
}