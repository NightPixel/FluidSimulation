#include "OpenGLUtils.h"
#include <vector>
#include <fstream>
#include <sstream>

// Reads all lines of text out of a file, and returns a string
// containing all those lines, separated by newline characters.
std::string readFile(const std::string& filePath)
{
    std::ifstream inStream(filePath);
    std::stringstream buffer;
    buffer << inStream.rdbuf();
    return buffer.str();
}

// Checks whether a shader compiled successfully.
// Returns (true, [info log string with warnings]) if compilation was successful.
// Otherwise, (false, [info log string with warnings/errors]) is returned.
std::pair<bool, std::string> checkShaderCompilation(GLuint shaderNumber)
{
    GLint status;
    glGetShaderiv(shaderNumber, GL_COMPILE_STATUS, &status);

    GLint logSize;
    glGetShaderiv(shaderNumber, GL_INFO_LOG_LENGTH, &logSize);

    std::vector<GLchar> logBuffer(logSize);
    glGetShaderInfoLog(shaderNumber, logSize, nullptr, logBuffer.data());

    if (status == GL_TRUE)
        return { true, logSize == 0 ? "" : logBuffer.data() };
    else
        return { false, logBuffer.data() };
}

// Creates a shader of the given type from the source code in the file at the given path.
GLuint createShaderFromSource(const std::string& sourceFilePath, GLenum shaderType)
{
    const std::string source = readFile(sourceFilePath);
    const char* sourceCharPtr = source.c_str();

    GLuint shader = glCreateShader(shaderType);
    glShaderSource(shader, 1, &sourceCharPtr, nullptr);
    glCompileShader(shader);

    return shader;
}
