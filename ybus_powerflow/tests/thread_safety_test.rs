// The types the C interface will hand out as handles are Send and Sync, as its thread-safety
//   contract (docs/FfiChecklist.md, phase F0) assumes: a handle can move between threads, and
//   one that isn't being changed can be read from many at once. The test passes by compiling.

use ybus_powerflow::ac::{AcModel, AcOptions, AcSolution};
use ybus_powerflow::dc::{DcModel, DcOptions, DcSolution, DcStudy};
use ybus_powerflow::opf::ac::{AcOpfOptions, AcOpfSolution};
use ybus_powerflow::opf::dc::{DcOpfOptions, DcOpfSolution};
use ybus_powerflow::{Case, PowerflowError};

fn assert_send_sync<T: Send + Sync>() {}

#[test]
fn the_types_behind_the_c_handles_are_send_and_sync() {
    assert_send_sync::<Case>();
    assert_send_sync::<DcModel>();
    assert_send_sync::<AcModel>();
    assert_send_sync::<DcStudy<'static>>();
    assert_send_sync::<DcOptions>();
    assert_send_sync::<AcOptions>();
    assert_send_sync::<DcOpfOptions>();
    assert_send_sync::<AcOpfOptions>();
    assert_send_sync::<DcSolution>();
    assert_send_sync::<AcSolution>();
    assert_send_sync::<DcOpfSolution>();
    assert_send_sync::<AcOpfSolution>();
    assert_send_sync::<PowerflowError>();
}
