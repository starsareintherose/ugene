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

#ifndef _U2_ADD_REFERENCE_TO_SEQUENCE_TASK_H_
#define _U2_ADD_REFERENCE_TO_SEQUENCE_TASK_H_

#include <U2Core/Task.h>
#include <U2Core/global.h>

#include <QMap>

namespace U2 {

class AssemblyObject;
class Document;
class GObjectSelection;
class LoadDocumentTask;
class U2SequenceObject;

/**
 * @AddReferencesToAssembliesByNameTask class finds reference sequence to one of selected assemblies.
 * Reference sequence should has the same name as the assembly scaffold has.
 */
class U2VIEW_EXPORT AddReferencesToAssembliesByNameTask : public Task {
public:
    AddReferencesToAssembliesByNameTask(const GObjectSelection* projectSelection, const QStringList& fileUrls);

    void prepare() override;
    QList<Task*> onSubTaskFinished(Task* subtask) override;
    ReportResult report() override;

private:
    bool findReferenceSequencesInDocument(Document* doc);
    void setAllReferences() const;

    const GObjectSelection* projectSelection = nullptr;
    QStringList fileUrls;

    QList<Task*> loadTasks;
    QList<Task*> addTasks;
    QList<Document*> loadedDocuments;
    QStringList notFounded;

    QMap<AssemblyObject*, U2SequenceObject*> assemblySequenceMap;
};

}


#endif