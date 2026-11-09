/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2026 team free-astro (see more in AUTHORS file)
 * Reference site is https://siril.org
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

#include <glib.h>
#include <glib/gstdio.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/mman.h>
#include <unistd.h>
#endif

#include "lsf/lsf_spill.h"

struct lsf_spill {
	void *data;
#ifdef _WIN32
	HANDLE file, map;
#else
	size_t bytes;
#endif
};

void *lsf_spill_data(const lsf_spill *sp) {
	return sp->data;
}

#ifdef _WIN32

/* removed by the system when the last handle closes, however Siril stops */
static HANDLE create_temp_file(const char *dir) {
	static gint serial = 0;
	HANDLE file = INVALID_HANDLE_VALUE;
	for (int attempt = 0; attempt < 100 && file == INVALID_HANDLE_VALUE; attempt++) {
		gchar *name = g_strdup_printf("siril_lsf_%lu_%d.tmp", GetCurrentProcessId(), g_atomic_int_add(&serial, 1));
		gchar *path = g_build_filename(dir, name, NULL);
		wchar_t *wpath = g_utf8_to_utf16(path, -1, NULL, NULL, NULL);
		if (wpath)
			file = CreateFileW(wpath, GENERIC_READ | GENERIC_WRITE, 0, NULL, CREATE_NEW,
					FILE_ATTRIBUTE_TEMPORARY | FILE_FLAG_DELETE_ON_CLOSE, NULL);
		gboolean again = file == INVALID_HANDLE_VALUE && GetLastError() == ERROR_FILE_EXISTS;
		g_free(wpath);
		g_free(path);
		g_free(name);
		if (file == INVALID_HANDLE_VALUE && !again)
			break;
	}
	return file;
}

lsf_spill *lsf_spill_new(const char *dir, size_t bytes) {
	lsf_spill *sp = g_new0(lsf_spill, 1);
	sp->file = create_temp_file(dir);
	LARGE_INTEGER size;
	size.QuadPart = (LONGLONG) bytes;
	if (sp->file != INVALID_HANDLE_VALUE
			&& (sp->map = CreateFileMappingW(sp->file, NULL, PAGE_READWRITE, size.HighPart, size.LowPart, NULL)))
		sp->data = MapViewOfFile(sp->map, FILE_MAP_ALL_ACCESS, 0, 0, bytes);
	if (!sp->data) {
		lsf_spill_free(sp);
		return NULL;
	}
	return sp;
}

void lsf_spill_free(lsf_spill *sp) {
	if (!sp)
		return;
	if (sp->data)
		UnmapViewOfFile(sp->data);
	if (sp->map)
		CloseHandle(sp->map);
	if (sp->file != INVALID_HANDLE_VALUE)
		CloseHandle(sp->file);
	g_free(sp);
}

#else

lsf_spill *lsf_spill_new(const char *dir, size_t bytes) {
	gchar *path = g_build_filename(dir, "siril_lsf_XXXXXX", NULL);
	int fd = g_mkstemp(path);
	if (fd < 0) {
		g_free(path);
		return NULL;
	}
	/* no name from here on: removed by the system when closed, however Siril stops */
	g_unlink(path);
	g_free(path);
	void *data = MAP_FAILED;
	if (!ftruncate(fd, (off_t) bytes))
		data = mmap(NULL, bytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	close(fd);	// the mapping keeps the file
	if (data == MAP_FAILED)
		return NULL;
	lsf_spill *sp = g_new0(lsf_spill, 1);
	sp->data = data;
	sp->bytes = bytes;
	return sp;
}

void lsf_spill_free(lsf_spill *sp) {
	if (!sp)
		return;
	munmap(sp->data, sp->bytes);
	g_free(sp);
}

#endif
