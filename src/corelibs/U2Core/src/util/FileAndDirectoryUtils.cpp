/**
 * UGENE - Integrated Bioinformatics Tools.
 * Copyright (C) 2008-2023 UniPro <ugene@unipro.ru>
 * http://ugene.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 */

#include "FileAndDirectoryUtils.h"

#include <QDir>
#include <QTemporaryFile>
#include <QTextStream>

#include <U2Core/GUrlUtils.h>
#include <U2Core/Log.h>
#include <U2Core/TextUtils.h>
#include <U2Core/TmpDirChecker.h>
#include <U2Core/U2SafePoints.h>

#ifdef Q_OS_WIN
#include <windows.h>

#endif

static const QString OUTPUT_SUBDIR = "run";

namespace U2 {

int FileAndDirectoryUtils::MIN_LENGTH_TO_WRITE = 32768;

#if defined(Q_OS_LINUX) | defined(Q_OS_DARWIN)
const QString FileAndDirectoryUtils::HOME_DIR_IDENTIFIER = "~/";
#else
const QString FileAndDirectoryUtils::HOME_DIR_IDENTIFIER = "%UserProfile%/";
#endif

QString FileAndDirectoryUtils::getFormatId(const FormatDetectionResult& r) {
    if (nullptr != r.format) {
        return r.format->getFormatId();
    }
    if (nullptr != r.importer) {
        return r.importer->getId();
    }
    return "";
}

QString FileAndDirectoryUtils::getWorkingDir(const QString& fileUrl, int dirMode, const QString& customDir, const QString& workingDir) {
    QString result;

    bool useInternal = false;

    if (dirMode == FILE_DIRECTORY) {
        result = GUrl(fileUrl).dirPath() + "/";
    } else if (dirMode == CUSTOM) {
        if (!customDir.isEmpty()) {
            result = customDir;
            if (!result.endsWith("/")) {
                result += "/";
            }
        } else {
            algoLog.error("Result folder is empty, default workflow folder is used");
            useInternal = true;
        }
    } else {
        useInternal = true;
    }

    if (useInternal) {
        result = workingDir;
        if (!result.endsWith("/")) {
            result += "/";
        }
        result += OUTPUT_SUBDIR;
        if (!result.endsWith("/")) {
            result += "/";
        }
        if (dirMode == WORKFLOW_INTERNAL_CUSTOM) {
            if (!customDir.isEmpty()) {
                result += customDir;
            }
            if (!result.endsWith("/")) {
                result += "/";
            }
        }
    }
    return result;
}

QString FileAndDirectoryUtils::createWorkingDir(const QString& fileUrl, int dirMode, const QString& customDir, const QString& workingDir) {
    const QString result = getWorkingDir(fileUrl, dirMode, customDir, workingDir);
    QDir dir(result);
    if (!dir.exists(result)) {
        dir.mkdir(result);
    }
    return result;
}

QString FileAndDirectoryUtils::detectFormat(const QString& url) {
    FormatDetectionConfig cfg;
    cfg.bestMatchesOnly = false;
    cfg.useImporters = true;
    cfg.excludeHiddenFormats = false;

    const QList<FormatDetectionResult> formats = DocumentUtils::detectFormat(url, cfg);
    if (formats.empty()) {
        return "";
    }

    return getFormatId(formats.first());
}

bool FileAndDirectoryUtils::isFileEmpty(const QString& url) {
    QFile file(url);
    if (!file.exists()) {
        return true;
    }
    if (file.size() == 0) {
        return true;
    }
    return false;
}

void FileAndDirectoryUtils::dumpStringToFile(QFile* f, QString& str) {
    CHECK(f != nullptr && str.length() >= MIN_LENGTH_TO_WRITE, );
    f->write(str.toLocal8Bit());
    str.clear();
}

bool FileAndDirectoryUtils::storeTextToFile(const QString& filePath, const QString& text) {
    QFile file(filePath);
    if (!file.open(QFile::WriteOnly)) {
        return false;
    }
    QTextStream out(&file);
    out << text;
    return true;
}

QString FileAndDirectoryUtils::getAbsolutePath(const QString& filePath) {
    CHECK(!filePath.isEmpty(), filePath);
    QString result = QDir::fromNativeSeparators(filePath);
    if (result.startsWith(HOME_DIR_IDENTIFIER, Qt::CaseInsensitive)) {
        result.remove(0, HOME_DIR_IDENTIFIER.length() - 1);
        result.prepend(QDir::homePath());
    }
    return QFileInfo(result).absoluteFilePath();
}

bool FileAndDirectoryUtils::isDirectoryWritable(const QString& dirPath) {
    QDir dir(dirPath);
    if (!dir.exists()) {
        return false;
    }

    QFile tmpFile(TmpDirChecker::getNewFilePath(dir.absolutePath(), "checkWritePermissions"));
    if (!tmpFile.open(QIODevice::WriteOnly)) {
        return false;
    }

    tmpFile.close();
    tmpFile.remove();
    return true;
}

bool FileAndDirectoryUtils::canWriteToPath(const QString& absoluteDirPath) {
    if (absoluteDirPath.isEmpty()) {
        return true;
    }
    QFileInfo fi(absoluteDirPath);
    SAFE_POINT(fi.dir().isAbsolute(), "Not an absolute path!", false);

    // Find out the folder that exists
    QDir existenDir(absoluteDirPath);
    while (!existenDir.exists()) {
        // Get upper folder
        QString dirPath = existenDir.path();
        QString dirName = existenDir.dirName();
        dirPath.remove(  // remove dir name and slash (if any) from the path
            dirPath.length() - dirName.length() - 1,
            dirName.length() + 1);
        if (dirPath.isEmpty()) {
            return false;
        }
        existenDir.setPath(dirPath);
    }

    // Attempts to write a file to the folder.
    // This assumes possibility to create any sub-folder, file, etc.
    QFile file(existenDir.filePath("testWriteAccess.txt"));
    if (!file.open(QIODevice::WriteOnly)) {
        return false;
    }
    file.close();
    file.remove();

    return true;
}

NP<FILE> FileAndDirectoryUtils::openFile(const QString& fileUrl, const QString& mode) {
#ifdef Q_OS_WIN
    QScopedPointer<wchar_t> unicodeFileName(TextUtils::toWideCharsArray(GUrlUtils::getNativeAbsolutePath(fileUrl)));
    QString modeWithBinaryFlag = mode;
    if (!modeWithBinaryFlag.contains("b")) {
        modeWithBinaryFlag += "b";  // Always open file in binary mode, so any kind of sam, sam.gz, bam, bai files are processed the same way.
    }
    QScopedPointer<wchar_t> unicodeMode(TextUtils::toWideCharsArray(modeWithBinaryFlag));
    auto err = errno;
    FILE* fp = _wfopen(unicodeFileName.data(), unicodeMode.data());
    err = errno;
    coreLog.error(QString("fopen_s: err: \"%1\"").arg(err));
    DWORD errorMessageID = ::GetLastError();
    if (errorMessageID != 0) {
        LPSTR messageBuffer = nullptr;

        //Ask Win32 to give us the string version of that message ID.
        //The parameters we pass in, tell Win32 to create the buffer that holds the message for us (because we don't yet know how long the message string will be).
        size_t size = FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
            NULL, errorMessageID, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPSTR)&messageBuffer, 0, NULL);

        //Copy the error message into a std::string.
        QString localErr = QString::fromUtf8(std::string(messageBuffer, size).c_str());
        coreLog.error(QString("fopen_s: can't open file: \"%1\", opening mode: %2, code: %3, error: %4, fp %5 NULL").arg(fileUrl).arg(mode).arg(errorMessageID).arg(localErr).arg(fp == nullptr ? "is" : "is not"));

        //Free the Win32's string's buffer.
        LocalFree(messageBuffer);
    }
    return fp;
#else
    return fopen(fileUrl.toLocal8Bit(), mode.toLatin1());
#endif
}

/** Closes file descriptor if the file descriptor is defined and is open. */
void FileAndDirectoryUtils::closeFileIfOpen(FILE* file) {
    int fd = file == nullptr ? -1 : fileno(file);
    CHECK(fd > 0, );
    fclose(file);
}

}  // namespace U2
