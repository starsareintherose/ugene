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

#include "AddReferencesToAssembliesByNameTask.h"

#include <U2Core/AddDocumentTask.h>
#include <U2Core/AssemblyObject.h>
#include <U2Core/AppContext.h>
#include <U2Core/DNASequenceObject.h>
#include <U2Core/DocumentModel.h>
#include <U2Core/GObject.h>
#include <U2Core/GObjectSelection.h>
#include <U2Core/GObjectTypes.h>
#include <U2Core/L10n.h>
#include <U2Core/LoadDocumentTask.h>
#include <U2Core/U2SafePoints.h>

#include <U2Gui/ProjectUtils.h>
#include <U2Gui/OpenViewTask.h>

#include <QList>

namespace U2 {

AddReferencesToAssembliesByNameTask::AddReferencesToAssembliesByNameTask(const GObjectSelection* _projectSelection,
                                                                         const QStringList& _fileUrls) :
    Task(tr("Add references to assemblies by name task"), TaskFlags_NR_FOSE_COSC),
    projectSelection(_projectSelection),
    fileUrls(_fileUrls) {}

void AddReferencesToAssembliesByNameTask::prepare() {
    const auto& selectedObjects = projectSelection->getSelectedObjects();
    for (auto obj : qAsConst(selectedObjects)) {
        auto assembly = qobject_cast<AssemblyObject*>(obj);
        SAFE_POINT(assembly != nullptr, L10N::nullPointerError("AssemblyObject"), );

        assemblySequenceMap.insert(assembly, nullptr);
    }

    for (const auto& filePath : qAsConst(fileUrls)) {
        if (ProjectUtils::hasLoadedDocument(filePath)) {
            auto document = ProjectUtils::findDocument(filePath);
            SAFE_POINT(document != nullptr, L10N::nullPointerError("Document"), );

            findReferenceSequencesInDocument(document);
        } else {
            auto loadDocTask = LoadDocumentTask::getDefaultLoadDocTask(stateInfo, filePath);
            CHECK_OP(stateInfo, );

            addSubTask(loadDocTask);
            loadTasks << loadDocTask;
        }
    }
    CHECK(loadTasks.isEmpty(), );

    setAllReferences();
}

QList<Task*> AddReferencesToAssembliesByNameTask::onSubTaskFinished(Task* subtask) {
    CHECK_OP(stateInfo, {});

    QList<Task*> res;
    if (loadTasks.contains(subtask)) {
        auto loadTask = qobject_cast<LoadDocumentTask*>(subtask);
        SAFE_POINT(loadTask != nullptr, L10N::nullPointerError("LoadDocumentTask"), {});

        auto loadedDocument = loadTask->takeDocument();
        SAFE_POINT(loadedDocument != nullptr, L10N::nullPointerError("Document"), {});

        loadedDocuments << loadedDocument;
        bool found = findReferenceSequencesInDocument(loadedDocument);
        if (!found) {
            loadedDocuments.removeOne(loadedDocument);
            loadedDocument->deleteLater();
        }
        loadTasks.removeOne(loadTask);
        CHECK(loadTasks.isEmpty(), {});

        for (auto ld : qAsConst(loadedDocuments)) {
            res << new AddDocumentTask(ld);
        }
        addTasks = res;
    } else if (addTasks.contains(subtask)) {
        addTasks.removeOne(subtask);
        CHECK(addTasks.isEmpty(), {});

        setAllReferences();
    }

    return res;
}

Task::ReportResult AddReferencesToAssembliesByNameTask::report() {
    CHECK_OP(stateInfo, Task::ReportResult::ReportResult_Finished);

    if (fileUrls.size() == notFounded.size()) {
        setError(tr("No of set files has sequence with a name, simular to selected assemblies names"));
    } else if (!notFounded.isEmpty()) {
        stateInfo.addWarning(tr("The following files do not contain sequences with a name, simular to "
            "selected assemblies names: %1").arg(notFounded.join(", ")));
    }

    return Task::ReportResult::ReportResult_Finished;
}

bool AddReferencesToAssembliesByNameTask::findReferenceSequencesInDocument(Document* doc) {
    auto docObjects = doc->getObjects();
    bool res = false;
    for (auto obj : qAsConst(docObjects)) {
        CHECK_CONTINUE(obj->getGObjectType() == GObjectTypes::SEQUENCE);

        auto casted = qobject_cast<U2SequenceObject*>(obj);
        SAFE_POINT(casted != nullptr, L10N::nullPointerError("U2SequenceObject"), false);

        auto assemblies = assemblySequenceMap.keys();
        for (auto assembly : qAsConst(assemblies)) {
            CHECK_CONTINUE(assemblySequenceMap.value(assembly) == nullptr);
            CHECK_CONTINUE(assembly->getGObjectName() == casted->getSequenceName());

            assemblySequenceMap.insert(assembly, casted);
            res = true;
            break;
        }
    }
    if (!res) {
        notFounded << doc->getURL().getURLString();
    }

    return res;
}

void AddReferencesToAssembliesByNameTask::setAllReferences() const {
    auto assemblies = assemblySequenceMap.keys();
    for (auto assembly : qAsConst(assemblies)) {
        auto referenceSequence = assemblySequenceMap.value(assembly);
        CHECK_CONTINUE(referenceSequence != nullptr);

        assembly->emitSetReference(referenceSequence);
    }
}



}