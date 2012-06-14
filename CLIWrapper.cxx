#ifdef WIN32
#define MODULE_IMPORT __declspec(dllimport)
#else
#define MODULE_IMPORT
#endif

extern "C" MODULE_IMPORT int ModuleEntryPoint(int, char* []);

int main(int argc, char** argv)
{
  return ModuleEntryPoint(argc, argv);
}
