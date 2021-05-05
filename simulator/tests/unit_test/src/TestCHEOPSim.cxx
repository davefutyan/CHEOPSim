#define BOOST_TEST_MODULE TestCHEOPSim
#include "boost/test/unit_test.hpp"
#include "boost/test/unit_test_suite.hpp"
#include "boost/test/unit_test_log.hpp"
#include "boost/test/test_tools.hpp"
#include "boost/test/detail/global_typedef.hpp"
#include "boost/filesystem.hpp"

#include <pqxx/pqxx>
#include <stdlib.h>

#define private public
#include "simulator/include/Simulator.hxx"
#include "ProgramParams.hxx"

using namespace boost;
using namespace boost::unit_test;

/// @brief unit test fixture
struct SimulatorFixture {
	SimulatorFixture() {
		m_params = CheopsInit(framework::master_test_suite().argc,
						      framework::master_test_suite().argv);

		m_simulator = new Simulator(m_params);
		m_data = m_simulator->m_data;

	}
	~SimulatorFixture() {
		string cleanup = "rm "+string(getenv("PWD"))+"/"+m_data->getOutputDirectory()+"/*/*.fits";
		std::system(cleanup.c_str());
		delete m_simulator;
	}
	ParamsPtr m_params; ///< ProgramParams pointer containing the input parameters
	Simulator * m_simulator; ///< Pointer to the simulator
	Data * m_data; ///< Pointer to the data
};

BOOST_FIXTURE_TEST_SUITE( TestCHEOPSim, SimulatorFixture )

BOOST_AUTO_TEST_CASE( testAll )
{

    m_simulator->process();

    BOOST_CHECK( m_data->stackedImageCount() ==  m_simulator->m_data->getTimeConfiguration().getNumberOfStackedImagesPITL());

}

BOOST_AUTO_TEST_SUITE_END()
