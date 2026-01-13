const { loadPyodide } = require("pyodide");
const { glob } = require("glob");

module.exports = async () => {
  let pyodide = await loadPyodide({
    indexURL: "node_modules/pyodide",
  });
  await pyodide.loadPackage("numpy")
  const wheels = await glob("dist/*whl");
  for (const wheel of wheels) {
    await pyodide.loadPackage(wheel); }
  console.log("Loaded Pyodide and wheels");
  pyodide.mountNodeFS("/mnt", "tests");
  global.pyodide = pyodide;
};
