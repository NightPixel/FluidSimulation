#pragma once

#include <string>
#include <GL/glew.h>

// Reads all lines of text out of a file, and returns a string
// containing all those lines, separated by newline characters.
std::string readFile(const std::string& filePath);

// Checks whether a shader compiled successfully.
// Returns (true, [info log string with warnings]) if compilation was successful.
// Otherwise, (false, [info log string with warnings/errors]) is returned.
std::pair<bool, std::string> checkShaderCompilation(GLuint shaderNumber);

// Creates a shader of the given type from the source code in the file at the given path.
GLuint createShaderFromSource(const std::string& sourceFilePath, GLenum shaderType);
