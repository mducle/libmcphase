beforeAll(() => {
  global.pyodide.runPython(`
    import importlib.util
    from importlib.machinery import SourceFileLoader
    import os
    import sys
    import unittest
    from unittest.mock import patch
    
    def run_test(testname):
        pathname = f"/mnt/{testname}.py"
        if not os.path.exists(pathname):
            raise ValueError(f"Test file {testname} not found")
        test_module_name = os.path.splitext(os.path.basename(pathname))[0]
        test_loader = SourceFileLoader(test_module_name, pathname)
        test_spec = importlib.util.spec_from_loader(test_module_name, test_loader)
        test_module = importlib.util.module_from_spec(test_spec)
        test_loader.exec_module(test_module)
        test_module_globals = dir(test_module)
        this_globals = globals()
        for key in test_module_globals:
            this_globals[key] = getattr(test_module, key)
    
        # create runner & execute
        print(f'Running {testname}')
        with patch('sys.exit') as exit_call:
            unittest.main(
                module=test_module,
                # We've processed the test source so don't let unittest try to reparse it
                # This forces it to load the tests from the supplied module
                argv=(sys.argv[0],),
                # these make sure that some options that are not applicable
                # remain hidden from the help menu.
                failfast=False,
                buffer=False,
                catchbreak=False,
            )
            # On success, unittest exits with code 0. Code 1 on failure
            exit_call.assert_called_with(0)
  `);
});

test("test_cfpars.py", () => { let result = global.pyodide.runPython(` run_test("test_cfpars") `); });
test("test_cf1ion.py", () => { let result = global.pyodide.runPython(` run_test("test_cf1ion") `); });
test("test_ic1ion.py", () => { let result = global.pyodide.runPython(` run_test("test_ic1ion") `); });
