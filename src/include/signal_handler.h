#ifndef _CONSOLE_HANDLER_H
#define _CONSOLE_HANDLER_H

#ifdef _WIN32
#include <Windows.h>
#else
#include <signal.h>
#endif

#include<functional>

namespace {

std::function<void(int)> _console_handler_ptr = nullptr;

#ifdef _WIN32
static BOOL WINAPI console_handler_internal(DWORD signal)
{
    if(_console_handler_ptr) {
        _console_handler_ptr(signal);
    }
    return TRUE;
}
#else
static void signal_handler_internal(int signal)
{
    if(_console_handler_ptr) {
        _console_handler_ptr(signal);
    }
}
#endif

bool set_signal_handler(std::function<void(int)> handler)
{
    _console_handler_ptr = handler;

#ifdef _WIN32
    DWORD console_flags = 0;
    HANDLE input_handle = GetStdHandle(STD_INPUT_HANDLE);
    GetConsoleMode(input_handle, &console_flags);
    SetConsoleMode(input_handle, ENABLE_EXTENDED_FLAGS | (console_flags & ~ENABLE_QUICK_EDIT_MODE));

    if(!SetConsoleCtrlHandler(console_handler, TRUE)) {
        return false;
    }
#else
    struct sigaction sig_handler;
    sig_handler.sa_handler = signal_handler_internal;
    sigemptyset(&sig_handler.sa_mask);
    sig_handler.sa_flags = 0;
    sigaction(SIGINT, &sig_handler, NULL);
    sigaction(SIGTERM, &sig_handler, NULL);

    // When mpirun receives SIGINT or SIGTERM, it forwards them to the MPI processes but then
    // exits before the MPI processes can exit cleanly. So instead use SIGUSR1. mpirun will forward
    // this to the MPI processes but will not exit early.
    sigaction(SIGUSR1, &sig_handler, NULL);
#endif

    return true;
}

};
#endif