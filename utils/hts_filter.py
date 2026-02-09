"""Suppress repeated htslib BAM-index warnings (show at most once in main process)."""
import sys
import os
import multiprocessing
import threading

_HTS_IDX_WARNING = "hts_idx_load3"
_OLDER_THAN_DATA = "older than the data file"
_warning_seen = False
_warning_lock = threading.Lock()


def _is_main_process():
    """Check if we're in the main process."""
    try:
        return multiprocessing.current_process().name == "MainProcess"
    except Exception:
        return False


def _filter_stderr_thread(pipe_r, real_stderr_fd):
    """Thread that reads from pipe and filters htslib warnings."""
    global _warning_seen
    real_stderr = os.fdopen(real_stderr_fd, "w", buffering=1)
    buffer = b""
    try:
        while True:
            try:
                data = os.read(pipe_r, 4096)
                if not data:
                    break
                buffer += data
                # Process complete lines
                while b"\n" in buffer:
                    line_bytes, buffer = buffer.split(b"\n", 1)
                    line = line_bytes.decode("utf-8", errors="replace") + "\n"
                    if _HTS_IDX_WARNING in line and _OLDER_THAN_DATA in line:
                        with _warning_lock:
                            if _is_main_process() and not _warning_seen:
                                _warning_seen = True
                                real_stderr.write(line)
                                real_stderr.flush()
                    else:
                        real_stderr.write(line)
                        real_stderr.flush()
            except OSError:
                break
        # Flush remaining buffer
        if buffer:
            line = buffer.decode("utf-8", errors="replace")
            if not (_HTS_IDX_WARNING in line and _OLDER_THAN_DATA in line):
                real_stderr.write(line)
                real_stderr.flush()
    except Exception:
        pass
    finally:
        try:
            real_stderr.close()
        except Exception:
            pass


def install_hts_warning_filter():
    """Install filter at file descriptor level to catch htslib warnings."""
    global _warning_seen
    if hasattr(install_hts_warning_filter, "_installed"):
        return
    install_hts_warning_filter._installed = True
    
    try:
        # Save original stderr file descriptor (fd 2)
        original_stderr_fd = sys.stderr.fileno()
        # Create a copy of the original stderr fd before redirecting
        saved_stderr_fd = os.dup(original_stderr_fd)
        # Create pipe
        pipe_r, pipe_w = os.pipe()
        # Redirect stderr (fd 2) to pipe write end
        os.dup2(pipe_w, original_stderr_fd)
        os.close(pipe_w)
        # Reopen sys.stderr to use the redirected fd
        sys.stderr.close()
        sys.stderr = os.fdopen(original_stderr_fd, "w", buffering=1)
        # Start filter thread that reads from pipe and writes to saved original stderr
        thread = threading.Thread(
            target=_filter_stderr_thread,
            args=(pipe_r, saved_stderr_fd),
            daemon=True
        )
        thread.start()
    except Exception:
        # If redirection fails, fall back to sys.stderr wrapper
        if not getattr(sys.stderr, "_is_hts_filter", False):
            sys.stderr = _HtsWarningFilter(sys.stderr)


class _HtsWarningFilter:
    """Fallback filter for sys.stderr (used if fd redirection fails)."""
    _is_hts_filter = True

    def __init__(self, stream):
        self._stream = stream
        self._seen_warning = False

    def write(self, buf):
        if _HTS_IDX_WARNING in buf and _OLDER_THAN_DATA in buf:
            if _is_main_process() and not self._seen_warning:
                self._seen_warning = True
                self._stream.write(buf)
            return
        self._stream.write(buf)

    def flush(self):
        self._stream.flush()

    def __getattr__(self, name):
        return getattr(self._stream, name)
